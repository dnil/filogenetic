#!/usr/bin/perl -w
#
# Gene Structure Stats 
# 
# Daniel Nilsson, 2010-04-19
#

use strict;

use Bio::Seq;
use Bio::SeqIO;

my $DEBUG = 0;
my $WARNING = 1;

sub usage { 
    print "Usage: gene_structure_stats.pl gbkfilename\n";
}

if ( @ARGV == 0 ) {
    usage;
    exit;
}

my $gbk_file = $ARGV[0];

my $gbk	= Bio::SeqIO->new('-file' => "<$gbk_file", -format=>'embl');

my $genes = 0;
my $transcripts = 0;

my $exons = 0;
my $introns = 0;

my %gene_transcripts;

my %transcript_introns;
my %transcript_exons;

my $sequences = 0;
my $exon_length = 0;
my $intron_length = 0;
my $gene_length = 0;
my $sequence_length = 0;

my @sequence_length;
my @intron_length;
my @exon_length;
my @gene_length;

# parse file and count features

while (	my $seq = $gbk->next_seq() ) {

    $DEBUG && print "DEBUG: Counting ",$seq->display_id,"\n";

    $sequence_length += $seq->length();
    $sequences++;
    push @sequence_length, $seq->length();

    my @features = $seq->get_SeqFeatures(); # just top level
    foreach my $feat ( @features ) {
	$DEBUG && print "DEBUG: feat ",$feat->primary_tag, ".\n";
        if($feat->primary_tag eq 'gene') {
	    # found gene
	    $genes++;
	    $gene_length += $feat->length;
	    push(@gene_length,$feat->length);

            if( $feat->has_tag("ID") ) {
                my ($gene_id) = ($feat->get_tag_values('ID'));
		$DEBUG && print "DEBUG: Inspecting $gene_id.\n";
            }
        }
	
	if($feat->primary_tag eq 'transcript') {
	    my $gene_id;
	    # parent = g1; ID =
            if( $feat->has_tag('Parent') ) {
		($gene_id) = ($feat->get_tag_values('Parent'));
            }
            if( $feat->has_tag("ID") ) {
                my ($transcript_id) = ($feat->get_tag_values('ID'));
		if(defined($gene_id)) {
		    $gene_transcripts{$gene_id}++;
		    $DEBUG && print "DEBUG: Inspecting transcript $transcript_id for $gene_id.\n";
		}
	    }
	    # found transcript
	    $transcripts++;
	}

        # exons are not explicitly given from AUGUSTUS. , note that the 5pUTR+first CDS or 3pUTR+lastCDS could be one exon..
        if($feat->primary_tag eq 'CDS' or $feat->primary_tag eq 'five_prime_UTR' or $feat->primary_tag eq 'three_prime_UTR' ) {
            # parent =g1.t1
            if( $feat->has_tag('Parent') ) {
                my ($transcript_id) = ($feat->get_tag_values('Parent'));

		if (defined($transcript_exons{$transcript_id})) {
		    $transcript_exons{$transcript_id}++;
		} else {
		    $transcript_exons{$transcript_id}=1;
		}

		push @exon_length,$feat->length;
            }

	    $exon_length += $feat->length;
	    $exons++;
        }

	if($feat->primary_tag eq 'intron') {
            if( $feat->has_tag('Parent') ) {
                my ($transcript_id) = ($feat->get_tag_values('Parent'));

		if (defined($transcript_introns{$transcript_id})) {
		    $transcript_introns{$transcript_id}++;
		} else {
		    $transcript_introns{$transcript_id}=1;
		}
		
		push(@intron_length,$feat->length);
            }

	    $intron_length+=$feat->length;
	    $introns++;
        }
    }
}

$DEBUG && print "DEBUG && Count features..\n";

# calculate 
my $n50val = N50(@sequence_length);
my ($sequence_mean_length, $sequence_median_length) = (mean_and_median(@sequence_length));
my ($gene_mean_length, $gene_median_length) = (mean_and_median(@gene_length));
my ($exon_mean_length, $exon_median_length) = (mean_and_median(@exon_length));
my ($intron_mean_length, $intron_median_length) = (mean_and_median(@intron_length));

#$exons = scalar( @exon_length );
#$introns = scalar( @intron_length );
#$genes = scalar( @gene_length );

my ($mean_introns_per_transcript, $median_introns_per_transcript) = (mean_and_median(values(%transcript_introns)));
my ($mean_exons_per_transcript, $median_exons_per_transcript) = (mean_and_median(values(%transcript_exons)));

# print

print "total sequence length ", sprintf("%.1f",$sequence_length/1000000)," Mb, average sequence length ", sprintf("%.1f",$sequence_length/$sequences)," nt\n";
print "total gene length ", sprintf("%.1f",$gene_length/1000000), " Mb, average gene length ", sprintf("%.1f",$gene_length/$genes)," nt\n";
print "$genes genes, $transcripts transcripts\n";
print "gene density ", sprintf("%.1f", $genes / $sequence_length * 1000000), " genes/Mb or ",sprintf("%.1f", $gene_length / $sequence_length *100) ,"% seq gene covered\n";

my $exon_estimate = $introns + 1*$genes;
print "~$exon_estimate exons ($exons CDSs and UTRs), $introns introns\n";

print "total exon length ", sprintf("%.1f",$exon_length/1000000)," Mb, total intron length ", sprintf("%.1f",$intron_length/1000000), " Mb \n";
print "average exon length ", sprintf("%.1f",$exon_length/$exon_estimate)," nt, average intron length ", sprintf("%.1f",$intron_length/$introns), " nt\n";
print "average number of introns/gene ", sprintf("%.1f",$introns/$genes),  ", average number of exons/gene ", sprintf("%.1f",$exon_estimate/$genes),"\n";
print "\n";
print "sequences\t$sequences:mean length\t$sequence_mean_length:\tmedian length\t$sequence_median_length:N50\t$n50val.\n";
print "genes\t$genes:mean length\t$gene_mean_length:\tmedian length\t$gene_median_length.\n";
print "exons\t$exons:mean length\t$exon_mean_length:\tmedian length\t$exon_median_length.\n";
print "introns\t$introns:mean length\t$intron_mean_length:\tmedian length\t$intron_median_length.\n";
print "\n";
print "introns per transcript:mean\t$mean_introns_per_transcript:\tmedian\t$median_introns_per_transcript.\n";
print "exons per transcript:mean\t$mean_exons_per_transcript:\tmedian\t$median_exons_per_transcript.\n";

sub mean_and_median {

    use POSIX; 

    my $sum=0; 
    my $n=0;
    my @val = @_;

    map { $sum += $_; } @val;
    $n = scalar(@val);

    my $median;
    if ($n%2==0){ 
	$median = ($val[floor($n/2)-1]+$val[floor($n/2)])/2;
    } else { 
	$median = $val[$n/2-1];
    } 

    return(sprintf("%.1f",$sum/$n), $median);
}

sub N50 {
    use POSIX; 
    my @val = sort {$b <=> $a} @_; # sort num desc.
    
    my $sum = 0;
    map { $sum += $_; } @val;
    
    my $half = $sum/2;

    my $cumsum= 0;

    my $current;
    while ( $cumsum < $half ) {
	$current = shift @val;
	$cumsum += $current;
    }
   
    return $current;
}

