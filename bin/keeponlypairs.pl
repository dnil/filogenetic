#!/usr/bin/perl -w
#
# (c) Daniel Nilsson, 2010
#

my $row;
my $ncount = 0;
my $bad = 0;
my @r;

my $badcount=0;
my $total=0;

my $expectedmate = 1 ;
my $pair =1 ;

while ( $row = <STDIN> ) {
    if ($ncount == 0) {
	push @r,$row;
	$ncount++;
#	print "DEBUG: ", substr($row, -2);

	if ( substr ($row, -2) ne "$expectedmate\n" ) {
	    $bad = 1;
	    if ( $expectedmate == 2 ) {
		$expectedmate = 1;
	    }	    
	} else {
	    #ok! found correct, bad=0, swap expected 
	    if ($expectedmate == 1) {
		$expectedmate =2;
	    } else {
		$expectedmate =1;
	    }
	}
    } elsif ($ncount == 1 || $ncount==2 ) {
	push @r,$row;
	$ncount++;
    } elsif ($ncount == 3) {
	$ncount =0;
	$total++;
	if ($bad == 0) {
	    if( $pair == 2) {
		print @r, $row;
		$pair=1;
		@r =();
	    } else {
		# first mate ok
		push @r,$row;
		$pair=2;
	    }
	} else {
	    $badcount++;
	    @r=();
	    $pair=1;
	}
	$bad=0;
    }

}
print STDERR "evicted $badcount/$total unpaired.\n";
