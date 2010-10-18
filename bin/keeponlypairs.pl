#!/usr/bin/perl -w
#
# (c) Daniel Nilsson, 2010
#

my $row;
my $ncount = 0;

my @r;

my $badcount=0;
my $total=0;

my $expectedmate = 1;
my $potentially_paired = 0;

my $matename = "";

while ( $row = <STDIN> ) {
    if ($ncount == 0) {
	$ncount++;
#	print "DEBUG: ", substr($row, -2);

	if ( substr ($row, -2) eq "1\n" ) {
	    @r=($row);
		$expectedmate = 2;
	    if ( $expectedmate == 1 ) {
		# ok, starting new record.
	    } else { # expectedmate ==2 , but got a 1 read
		# replace the previously stored /1 record with the current, and set expectedmate to 2
	    }
	} elsif (substr ($row, -2) eq "2\n" ) {
	    push @r, $row;
	    if ($expectedmate == 1) {
                # found a /2 read when waiting for a /1
                # don't really care about @row,as this entry will be discarded.
		$badcount++;
	    } else { #expectedmate==2
		$expectedmate =1;
		$potentially_paired=1; #ok! found a potential pair! push row to memory for printing if it holds.
	    }
	    
	} else {
	    print "something is wrong with the read naming: got a PE named $row";
	    exit;
	}
    } elsif ($ncount == 1 || $ncount==2 ) {
	push @r,$row;
	$ncount++;
    } elsif ($ncount == 3) {
	$ncount =0;
	$total++;
	if( $potentially_paired == 1) {
	    
	    # check if really named ok
	    if ($matename eq substr($r[4],0,-2)) {
		# yes, so print.
		print @r, $row;
	    } else {
		# actually bad, although not yet marked so
		$badcount++;
	    }
	    #lower flag
	    $potentially_paired=0;
	    # clear entry..
	    @r =();
	} else { #potentially paired ==0
	    # first mate ok, proceed
	    push @r,$row;
#	    $potentially_paired=1; will be set by expectedmate check at row 0..
	    $matename = substr($r[0],0,-2);
	}
    }

}
print STDERR "evicted $badcount/$total unpaired.\n";
