#!/usr/bin/perl

=head1 NAME

contig_stats.pl - Takes one or more fasta files as input and calculates some assembly stats for it

=head1 SYNOPSIS

contig_stats.pl -fasta /path/to/file1.fasta -fasta /path/to/file2.fasta -cutoff 100 -cutoff 1000 -cutoff 10000 -delim ","

=head1 DESCRIPTION

contig_stats.pl takes one or more fasta files as input and calculates contig stats such as number, total assembled bases, N50, for each contig size cutoff. If you have R installed, it also gives you pretty plots (not yet)

=head1 AUTHORS

sujai.kumar@ed.ac.uk 2010.05.12

=cut

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);

my @fastafiles;
my @covfiles;
my @thresholds;
my $delim = "\t";
my $output_dir = "pc";
my $graphs = "";
my $length_cutoff;

GetOptions (
    "fastafile=s" => \@fastafiles,
    "covfile=s"   => \@covfiles,
    "threshold=i" => \@thresholds,
    "delimiter=s" => \$delim,
    "output=s"    => \$output_dir,
    "graphs"      => \$graphs,
    "length=i"    => \$length_cutoff,
);
@thresholds = (0,100,1000,10000) unless @thresholds;

#---------------------------------

die "Usage: contig_stats.pl -f contigs1.fa -c contigs1.fa.cov -f contigs2.fa -c contigs2.fa.cov -t 100 -t 1000\n-t threshold values are optional. -c cov files (corresponding to -f fasta files) are also optional and are tab delimited files where col 1 is contig id and col 2 is coverage\n" unless @fastafiles;

if (-e $output_dir) {
    print STDERR "$output_dir exists, will be overwritten\n";
    system "rm -rf $output_dir";
}
mkdir $output_dir or die "Could not create $output_dir\nCheck if you have permission\n";

my $toprint = "";

# print header of table
$toprint .= "Filename${delim}Max contig length";
foreach (@thresholds) {
    $toprint .= "${delim}Num contigs >$_${delim}Total bases in contigs >$_${delim}N50 for contigs >$_${delim}Contigs >${_} in N50${delim}GC contigs >$_${delim}nonATGC in contigs >$_${delim}Mean length for contigs >$_";
}
$toprint .= "\n";

my %sequences_all;
for my $fastafile (@fastafiles) {

    # Read in contigs from fasta files 
    $sequences_all{$fastafile} = &fastafile2hash($fastafile);
    my @sorted_contigs = 
        sort { $sequences_all{$fastafile}{$b}{len} <=> $sequences_all{$fastafile}{$a}{len} }
            keys % { $sequences_all{$fastafile} };
    my $longest_contig = $sequences_all{$fastafile}{$sorted_contigs[0]}{len};
    
    my $covfile;
    if (@covfiles) {
        $covfile = shift @covfiles;
        $sequences_all{$fastafile} = &addCoverage2Hash ($sequences_all{$fastafile}, $covfile);
    }
    &printContigLengths($output_dir, $fastafile, $sequences_all{$fastafile});
    
    $toprint .= "$fastafile${delim}$longest_contig";

    # for each cutoff 
    for my $threshold (@thresholds) {
        my $num_contigs = 0;
        my $total_bases = 0;
        my $N50 = 0;
        my $contigs_in_N50 = 0;
        my $gc_count = 0;
        my $total_nonatgc = 0;
        
        # calculate num contigs, total bases and gc count at this cutoff
        for my $contig (@sorted_contigs) {
            my $contig_len = $sequences_all{$fastafile}{$contig}{len};
            last if $contig_len < $threshold;
            $num_contigs++;
            $total_bases   += $contig_len;
            $total_nonatgc += $sequences_all{$fastafile}{$contig}{nonatgc};
            $gc_count      += $sequences_all{$fastafile}{$contig}{gc};
        }
        my $mean_contig_length;
        $mean_contig_length = $total_bases/$num_contigs if $num_contigs >0;

        # calculate N50
        my $cumulative_total = 0;
        for my $contig (@sorted_contigs) {
            my $contig_len = $sequences_all{$fastafile}{$contig}{len};
            $cumulative_total += $contig_len;
            $N50 = $contig_len;
            $contigs_in_N50 ++;
            last if ($cumulative_total > $total_bases/2);
        }
        
        if ($total_bases > 0) {
            $toprint .= "${delim}$num_contigs${delim}$total_bases${delim}$N50${delim}$contigs_in_N50${delim}" . 
                        sprintf("%.1f",$gc_count*100/($total_bases - $total_nonatgc + 1)) . 
                        ${delim} . $total_nonatgc . $delim . sprintf("%.1f",$mean_contig_length);
        }
        else {
            $toprint .= "${delim}0${delim}0${delim}NA${delim}NA${delim}NA${delim}NA";
        }
    }
    # new line at end of contig stats row for each fasta file
    $toprint .= "\n";
}
open  STAT, ">$output_dir/contig_stats.txt" or die $!;
print STAT $toprint;
print $toprint;
close STAT;

if ($graphs) {
open  TMP,">$output_dir/contig_stats.R" or die $!;
print TMP <<R;
pdf("$output_dir/contig_lengths_gc.pdf",10,10)
contigs=read.delim("$output_dir/contig_lengths_gc.txt",header=FALSE,col.names=c('file','id','length','gc','n','cov'))
files <- as.vector(unique(contigs\$file))
maxX = maxY = 0;
for (i in files) {
    s<-subset(contigs,file==i & length>=$graphs);
    l<-dim(s)[[1]];
    if (l > maxX) maxX = l;
    c<-sum(s\$length);
    if (c > maxY) maxY = c;
}
plot(0,xlim=c(0,maxX),ylim=c(0,maxY),cex=0.0001,xlab="Contigs ranked by size",ylab="Cumulative contig length")
j <- 0;
for (i in as.vector(unique(contigs\$file))) {
    j <- j+1;
    s <- subset(contigs,file==i & length>=$graphs)
    points(cumsum(sort(s\$length,decreasing=TRUE)),pch=".",col=j,cex=.5)
}
legend(0.5 * maxX, 0.5 * maxY, files, col=1:j,lty=c(1),lwd=3)
dev.off()
R
close TMP;
# run R script
system("R -f $output_dir/contig_stats.R >/dev/null 2>&1");
}

#############################################################################

sub fastafile2hash {
    my $fastafile = shift @_;
    my %sequences_lengc;
    my $fh = &read_fh($fastafile);
    my $header;
    while (my $line = <$fh>) {
        if ($line =~ /^>(\S+)(.*)/) {
            $header = $1;
            # $sequences{$header}{desc} = $2;
        }
        else {
            chomp $line;
            $sequences_lengc{$header}{len} += length $line;
            $sequences_lengc{$header}{gc}  += ($line =~ tr/gcGC/gcGC/);
            $line =~ s/[^atgc]/N/ig;
            $sequences_lengc{$header}{nonatgc} += ($line =~ tr/N/N/);
        }
    }
    close $fh;
    if ($length_cutoff) {
        foreach (keys %sequences_lengc) {
            delete $sequences_lengc{$_} if $sequences_lengc{$_}{len} < $length_cutoff;
        }
    }
    return \%sequences_lengc;
}

#############################################################################

sub addCoverage2Hash {
    my ($sequences_cov, $covfile) = @_;

    # takes coverage info for each contig from covfile, and attaches as cov field to each contig id key
    my $fh = &read_fh($covfile);
    while (my $line = <$fh>) {
        next unless $line =~ /^(\S+)\t(\S+)$/;
        $$sequences_cov{$1}{cov} = $2 if (defined $$sequences_cov{$1})
    }
    return $sequences_cov;
}

#############################################################################

sub printContigLengths {
    my ($output_dir, $fastafile, $sequences) = @_;
    
    open LEN, ">>$output_dir/contig_lengths_gc.txt" or die $!;
    foreach (keys %{$sequences}) {
        print LEN "$fastafile\t$_\t$$sequences{$_}{len}\t$$sequences{$_}{gc}\t$$sequences{$_}{nonatgc}";
        print LEN "\t$$sequences{$_}{cov}" if exists $$sequences{$_}{cov};
        print LEN "\n";
    }
    close LEN;
}

#############################################################################

sub read_fh {
    my $filename = shift @_;
    my $filehandle;
    if ($filename =~ /gz$/) {
        open $filehandle, "gunzip -dc $filename |" or die $!;
    }
    else {
        open $filehandle, "<$filename" or die $!;
    }
    return $filehandle;
}
