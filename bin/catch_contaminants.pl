#!/usr/bin/perl  -w 

my $Evalcutoff = 1e-5;
my $match_cutoff = 0.4;

my $contig_length;
my $contig;
my @matched;

my $last_contig="";
my $last_contig_length=0;
my $first=1;

my $DEBUG=0;

while(my $r=<STDIN>) {
    chomp $r;

    my @r = split (/\s+/,$r);
    
    $contig=$r[0];
    ($contig_length) = ($contig =~ m/countedlength\_(\d+)/);
    
    $DEBUG && print STDERR "$contig\t$contig_length\n";

    if( $contig eq $last_contig) {
	# same contig as last time around
	# store any high quality matches
	if ( $r[10] < $Evalcutoff ) {
	    
	    $start=$r[6] -1;
	    $end=$r[7] -1;
	    
	    $DEBUG && print STDERR "marking $start -- $end as matches since eval ",$r[10],"\n";

	    for (my $i = $start; $i <= $end ; $i++ ) {
		$matched[$i] = 1;
	    }
	}
    } else {
	# evaluate current contig
	if ($first == 0) {
	    my $nr_matched=0;
	    map { $nr_matched += $_ } @matched;
	
	    if ( $nr_matched / $last_contig_length > $match_cutoff ) {
		print "$last_contig\tCONTAMINANT\t$nr_matched\t$last_contig_length\n";
	    }
	}

	# next contig!
	@matched = ();
	@matched = (0)x$contig_length;
	$DEBUG && print STDERR "setting contig contig length ",scalar(@matched), " - should be $contig_length.\n";
	
	if ( $r[10] < $Evalcutoff ) {
	    
	    $start=$r[6] -1;
	    $end=$r[7] -1;
	    
	    $DEBUG && print STDERR "marking $start -- $end as matches since eval ",$r[10],"\n";

	    for (my $i = $start; $i <= $end ; $i++ ) {
		$matched[$i] = 1;
	    }
	}

	$last_contig=$contig;
	$last_contig_length=$contig_length;
	$first=0;
    }
}

# last contig..
my $nr_matched=0;
map { $nr_matched += $_ } @matched;

if ( $nr_matched / $contig_length > $match_cutoff ) {
    print "$contig\tCONTAMINANT\t$nr_matched\t$contig_length\n";
}
