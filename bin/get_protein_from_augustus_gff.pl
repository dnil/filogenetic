#!/usr/bin/perl -w 

my $in_seq =0;
my $gene_name = "missing";

while( my $r = <STDIN> ) {
    chomp $r;

    if($in_seq) {
	if( $r=~/^\#\s+(\w*)\]/) {
	    $seq .= $1;
	    
	    fasta_print_seq( $gene_name, $seq );

	    $in_seq = 0; 
	} elsif ($r=~/^\#\s+(\w+)/) {
	    $seq .= $1;
	}
	

    }
    
    if($r=~/\#\s+protein sequence \= \[(\w+)(\]?)/) {
	$in_seq = 1;
	$seq = $1;

	if($2 eq "]")
	{	    
	    fasta_print_seq( $gene_name, $seq );
	    $in_seq = 0;
	}
    }
# Contig_11_length_4832   AUGUSTUS        gene    541     4455    0.06    -       .       g8
    if($r =~ /^([^#]+)\s+AUGUSTUS\s+gene\s+(\d+)\s+(\d+)\s+([\d\.]+)\s+([-+]{1})/) {
	$gene_name = "$1_$2$5$3_$4";
    }
}


sub fasta_print_seq {
    my $gene_name =shift;
    my $seq = shift;

    print ">".$gene_name."\n";
    # prettyprint
    my @seq = split(/ */, $seq);		    
    while ( my @curr = splice @seq, 0, (80>=scalar(@curr))?80:scalar(@curr)) {
	print join("",@curr),"\n";
    }
    
}
