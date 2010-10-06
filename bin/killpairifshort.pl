#!/usr/bin/perl -w
#
# (c) Daniel Nilsson, 2010
#

my $row;
my $ncount = 0;
my $bad = 0;
my @r;

my $cutoff=35;

my $badcount=0;
my $total=0;

while ( $row = <STDIN> ) {
    if ($ncount == 0 || $ncount==3 ||$ncount==4) {
	$r[$ncount]=$row;
	$ncount++;
    } elsif ($ncount == 1 || $ncount==5) {
	if (length($row) < $cutoff+1) {
	    $bad=1;
	} else {
	    $r[$ncount] = $row;
	}
	$ncount++;
    } elsif ($ncount == 2 || $ncount==6) {
	$r[$ncount]=$row;
	$ncount++;
    } elsif ($ncount == 7) {
	$ncount = 0;

#	if ($row =~ m/B/) {
#	    $bad=1;
#	}

	$total++;
	if ($bad == 0) {
	    print @r, $row;
	} else {
	    $badcount++;
	}
	$bad=0;
    }

}
print STDERR "bad $badcount/$total pairs.\n";
