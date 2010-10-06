#!/usr/bin/perl -w
#
# (c) Daniel Nilsson, 2010
#

my $row;
my $ncount = 0;
my $bad = 0;
my $r0;
my $r1;
my $r2;
my $badcount=0;
my $total=0;

while ( $row = <STDIN> ) {
    if ($ncount == 0) {
	$r0=$row;
	$ncount++;
    } elsif ($ncount == 1) {
	if ($row =~m/N/) {
	    $bad=1;
	} else {
	    $r1 = $row;
	}
	$ncount++;
    } elsif ($ncount == 2) {
	$r2=$row;
	$ncount++;
    } elsif ($ncount == 3) {
	$ncount = 0;

#	if ($row =~ m/B/) {
#	    $bad=1;
#	}

	$total++;
	if ($bad ==0) {
	    print $r0, $r1, $r2, $row;
	} else {
	    $badcount++;
#	    print $r0, "N\n", $r2, "A\n";
	}
	$bad=0;
    }

}
print STDERR "bad $badcount/$total.\n";
