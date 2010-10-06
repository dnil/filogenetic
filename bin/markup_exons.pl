#!/usr/bin/perl -w
#
# (c) Daniel Nilsson, 2010-04-20
#
#
# Usage: Usage: markup_exons.pl -i gbkfile [-y output_file_type($def_output_file_type)] -o outputgenbankfile [-t output_file_type($def_output_file_type)] [-W|--warn] [-D|--debug]
#

=head1 NAME

markup_exons.pl

=head1 AUTHOR

Daniel Nilsson, daniel.nilsson@scilifelab.se, daniel.nilsson@ki.se, daniel.k.nilsson@gmail.com

=head1 LICENSE AND COPYRIGHT

Copyright 2009, 2010 held by Daniel Nilsson. The package is released for use under the Perl Artistic License.

=head1 SYNOPSIS

Usage: markup_exons.pl -i gbkfile [-y output_file_type($def_output_file_type)] -o outputgenbankfile [-t output_file_type($def_output_file_type)] [-W|--warn] [-D|--debug]

=cut

use strict;

use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;

my $DEBUG = 0;
my $WARNING = 1;

my $def_input_file_type = 'genbank';
my $def_output_file_type = 'genbank';

sub usage {
    print "Usage: markup_exons.pl -i gbkfile [-y output_file_type($def_output_file_type)] -o outputgenbankfile [-t output_file_type($def_output_file_type)] [-W|--warn] [-D|--debug]\n";
}

my $input_file_name = undef;
my $input_file_type = $def_input_file_type;
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

my $sin  = Bio::SeqIO->new(-file => "<$fasta_file_name", -format => $input_file_type);

my $sout = Bio::SeqIO->new(-file => ">$output_file_name" , -format => $output_file_type);

while ( my $seq = $sin->next_seq() ) {

    $DEBUG && print "DEBUG: Reading ",$seq->display_id,"\n";
    my $current_gene;
    my $current_transcript;

    foreach my $feat ( @features ) {
	$DEBUG && print "DEBUG: feat ",$feat->primary_tag, ".\n";
        if($feat->primary_tag eq 'gene') {
	    # found gene

            if( $feat->has_tag("ID") ) {
                 ($current_gene) = ($feat->get_tag_values('ID'));
		$DEBUG && print "DEBUG: Inspecting $gene_id.\n";

		push @genes, $gene_id;
            }

	    $current_transcript = undef;
        }

	if($feat->primary_tag eq 'transcript') {

	    # parent = g1; ID =
            if( $feat->has_tag('Parent') ) {
		($transcript_gene_id) = ($feat->get_tag_values('Parent'));
		if ($transcript_gene_id ne $current_gene_id) {
		    $WARNING && print "WARNING: found transcript with parent $transcript_gene_id in gene $current_gene_id. Unsorted input features?\n";
		}
            }
            if( $feat->has_tag("ID") ) {
                my ($transcript_id) = ($feat->get_tag_values('ID'));
		if(defined($gene_id)) {
		    $DEBUG && print "DEBUG: Inspecting transcript $transcript_id for $gene_id.\n";
		}
	    }
	    # found transcript
	    $transcripts++;
	}



}
