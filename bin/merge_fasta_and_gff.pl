#!/usr/bin/perl -w
#
# Daniel Nilsson, 2010-04-16
#

use strict;

use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;

my $DEBUG = 0;
my $WARNING = 1;

my $def_output_file_type = 'genbank';

sub usage {
    print "Usage: merge_fasta_and_gff.pl -g gfffile -f fastafile -o outputgenbankfile [-t output_file_type($def_output_file_type)] [-W|--warn] [-D|--debug]\n";
}

my $gff_file_name = undef;
my $fasta_file_name = undef;
my $output_file_name = undef;
my $output_file_type = $def_output_file_type;

while (my $arg = shift @ARGV) {
    if ($arg =~ /^-/) {
        if ($arg eq '-g') {
            my $next_arg = shift @ARGV;
            if($next_arg eq "") {
                print "-g requires an argument, but none given. Bailing out.\n";
		usage();
                exit 1;
            } else {
                $gff_file_name = $next_arg;    # A GENE GFF
            }
        }

        if ($arg eq '-f') {
            my $next_arg = shift @ARGV;
            if($next_arg eq "") {
                print "-f requires an argument, but none given. Bailing out.\n";
		usage();
                exit 1;
            } else {
                $fasta_file_name = $next_arg;    # A GENE GFF
            }
        }

        if ($arg eq '-o') {
            my $next_arg = shift @ARGV;
            if($next_arg eq "") {
                print "-o requires an argument, but none given. Bailing out.\n";
		usage();
                exit 1;
            } else {
                $output_file_name = $next_arg;
            }
        }

        if ($arg eq '-t') {
            my $next_arg = shift @ARGV;
            if($next_arg eq "") {
                print "-t requires an argument, but none given. Bailing out.\n";
		usage();
                exit 1;
            } else {
                $output_file_type = $next_arg;
            }
        }

	if ($arg eq '-W' or $arg eq '--warn') {
	    $WARNING = 1;
	}

	if ($arg eq '-D' or $arg eq '--debug') {
	    $DEBUG = 1;
	}
    }
}

if ( !defined($gff_file_name) or !defined($fasta_file_name) or !defined($output_file_name)) {
    print "Error. Required file name missing.\n";
    usage();
    exit 1;
}

# both fasta and gff are reasonably small footprint as text

my $sin  = Bio::SeqIO->new(-file => "<$fasta_file_name" ,
		       -format => 'Fasta');

my $sout = Bio::SeqIO->new(-file => ">$output_file_name" ,
		       -format => $output_file_type);

# hash with refrences to all seqs, indexed by primary_ID
my %seqid;

while ( my $seq = $sin->next_seq() ) {
    my $current_id = lc $seq->primary_id;

    if(defined($seqid{$current_id})) {
	$WARNING && print "WARNING: multiple entries in $fasta_file_name with id ", $seq->primary_id, ".\n";
    } else {
	# add an Annotation object to get seq name and such into gbk/embl-file?? Well, setting -accession_number or -id should be enough?!
	$seq->id( $current_id );
	$seq->accession_number( $current_id );

	$seq->display_id( $current_id ); # setting display_id is probably overkill??
	
	$seqid{$current_id}=\$seq;
	$DEBUG && print "Adding sequence $current_id ",$seq->accession_number(),".\n";
    }
}

open GFF, "<$gff_file_name" or die "Could not open GFF file $gff_file_name.\n";

while ( my $gff_line = <GFF> ) {
    
    my ($current_seq_id) = ($gff_line =~ m/^(\S+)\s+/);
    $current_seq_id = lc $current_seq_id;

    if ( $gff_line =~ m/^\s*\#/ ) {
	# ignore comment
    } elsif ( exists( $seqid{$current_seq_id} ) && defined( $seqid{$current_seq_id} ) ) {     
	
	my $feat = new Bio::SeqFeature::Generic();
	$feat->gff_format(Bio::Tools::GFF->new('-gff_version' => 3));
	$feat->set_attributes('-gff_string' => $gff_line);

	$DEBUG && print "Have seq matching feat line: $gff_line\n";
	$DEBUG && print "Adding feat ",$feat->primary_tag," to sequence $current_seq_id, with tags", join(" ",$feat->get_all_tags(),"\n");

	${$seqid{$current_seq_id}}->add_SeqFeature($feat);
    }
}

foreach my $seq ( keys %seqid ) {
    $DEBUG && print "Writing sequence ",${$seqid{$seq}}->accession_number(),"\n";
    $sout->write_seq(${$seqid{$seq}});
}
