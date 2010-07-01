#!/usr/bin/perl -w
# 
# Daniel Nilsson

my $inseq = 0; 
my $seq = "";
my $name = "";

while ($l = <STDIN>) {
    chomp $l;
    
    if ($inseq) {

	if ($l =~ m/^>(\S+)/) { 
	    $inseq=1; 

	    my $len = length($seq);
	    print ">".$name."_countedlength_".$len."\n";
	    # prettyprint	    
	    my @seq = split(/ */, $seq);		    
	    while ( my @curr = splice @seq, 0, (80>=scalar(@curr))?80:scalar(@curr)) {
		print join("",@curr),"\n";
	    }    	 

	    $name = $1;
	    $seq = "";
	} else {
	    $l =~ s/\s+//g;
	    $seq .= $l;	
	}	
    } elsif ($l =~ m/^>(\S+)/) { 
	# first time around..
	$inseq=1; 
	
	$name = $1;
	$seq = "";
    }
}

if($inseq==1) { 
    my $len = length($seq);
    print ">".$name."_".$len."_length_".$len."\n";
    # prettyprint

    my @seq = split(/ */, $seq);		    
    while ( my @curr = splice @seq, 0, (80>=scalar(@curr))?80:scalar(@curr)) {
	print join("",@curr),"\n";
    }    
}
