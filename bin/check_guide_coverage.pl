#!/usr/bin/perl  -w 

my $Evalcutoff = 1e-5;
my $match_cutoff = 0.8;
my $novel_bases_for_segment_cutoff = 0.6;

my $contig_length=0;
my $contig="";
my $hit="";
my @matched;
my $segments=0;

my $last_contig;
my $last_hit;
my $last_contig_length;
my $first=1;

my $DEBUG=0;

while(my $r=<STDIN>) {
    chomp $r;

    my @r = split (/\s+/,$r);

    $last_contig=$contig;
    $contig=$r[0];

    $last_hit=$hit;
    $hit=$r[1];

    $last_contig_length=$contig_length;
    ($contig_length) = ($contig =~ m/length\_(\d+)/);
    
    $DEBUG && print STDERR "$contig\t$contig_length\n";

    if( $contig eq $last_contig) {
	# same contig as last time around
	# store any high quality matches
	if ( $r[10] < $Evalcutoff ) {
	    
	    $start=$r[6] -1;
	    $end=$r[7] -1;
	    
	    $DEBUG && print STDERR "marking $start -- $end as matches since eval ",$r[10],"\n";

	    my $novel_bases =0;
	    for (my $i = $start; $i <= $end ; $i++ ) {	       
		if ( $matched[$i] == 0 ) {
		    $novel_bases++;
		}
		$matched[$i] = 1; 
	    }
	    if( $novel_bases / ($end-$start+1) > $novel_bases_for_segment_cutoff) {
		if( $hit ne $last_hit ) {
		    # disregard segmenting within a hit contig; e g exons for an EST..
		    # assumes hits are sorted, which they generally are.

		    $segments++;
		    $DEBUG && print STDERR "increasing to $segments segments since $hit ne $last_hit.\n";
		}
	    }
	}
    } else {
	# evaluate current (i e last) contig
	if ($first == 0) {
	    my $nr_matched=0;
	    map { $nr_matched += $_ } @matched;
	
	    if ( $nr_matched / $last_contig_length > $match_cutoff ) {
		print "$last_contig\tFOUND_OK\t$nr_matched\t$last_contig_length\t$segments\n";
	    } elsif( $nr_matched >0 ) {
		print "$last_contig\tHITS_ONLY\t$nr_matched\t$last_contig_length\t$segments\n";
	    } else {
		print "$last_contig\tNOT_FOUND\t$nr_matched\t$last_contig_length\t$segments\n";
	    }

	}

	# next contig!
	$segments=0;
	@matched = ();
	@matched = (0)x$contig_length;
	$DEBUG && print STDERR "setting contig contig length ",scalar(@matched), " - should be $contig_length.\n";
	
	if ( $r[10] < $Evalcutoff ) {
	    
	    $start=$r[6] -1;
	    $end=$r[7] -1;
	    
	    $DEBUG && print STDERR "marking $start -- $end as matches since eval ",$r[10],"\n";

	    my $novel_bases =0;
	    for (my $i = $start; $i <= $end ; $i++ ) {	
		if ( $matched[$i] == 0 ) {
		    $novel_bases++;
		}
		$matched[$i] = 1; 
	    }
	    if( $novel_bases / ($end-$start+1) > $novel_bases_for_segment_cutoff) {
	            # we dont check last hit name for the first hit for a contig, as that would be misleading or unneccesary.
		    # disregard segmenting within a hit contig; e g exons for an EST..
		    # assumes hits are sorted, which they generally are.
		    $segments++;
		    $DEBUG && print STDERR "increasing to $segments segments since $hit ne $last_hit.\n";
		    
	    }
	}

	$first=0;
    }
}

# last contig..
my $nr_matched=0;
map { $nr_matched += $_ } @matched;

if ( $nr_matched / $contig_length > $match_cutoff ) {
    print "$last_contig\tFOUND_OK\t$nr_matched\t$last_contig_length\t$segments\n";	       

} elsif( $nr_matched >0 ) {
    print "$last_contig\tHITS_ONLY\t$nr_matched\t$last_contig_length\t$segments\n";
} else {
    print "$last_contig\tNOT_FOUNDt\$nr_matched\t$last_contig_length\t$segments\n";
}
