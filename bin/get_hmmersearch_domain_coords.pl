#!/usr/bin/perl -w

my $indom= 0;

while(my $l = <STDIN>) {
    chomp $l;

    if($indom == 1) {
	if ($l=~/Sequence/ and $l=~/Domain/ and $l=~/seq-f/) {
	}
	elsif ($l=~/^--------\s+/ ) {
	}
	elsif ($l eq "") {
	    $indom = 0;
	    exit;
	}
	else {
	    @c=split(/\s+/,$l);
	    print $c[0],"\t", $c[2],"\t", $c[3], "\n";	    
	}       
    } elsif ($l =~ /^Parsed for domains:/) { 
	$indom=1;
    }
}

