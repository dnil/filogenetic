#!/usr/bin/perl  -w 

# Total number of bases and length
open (FILEHANDLE,"H:\DATA\PhD\Fasteris\Data\2009-12-07_GEU-1_velvet-h31_SCF.fasta");
my $contig=<FILEHANDLE>;
	my $A=($contig=~tr/A//);
	my $T=($contig=~tr/T//);
	my $G=($contig=~tr/G//);
	my $C=($contig=~tr/C//);
print"A:$A";
print"T:$T";
print"G:$G";
print"C:$C";

	my $DNA_num_bases=$A+$T+$G+$C;
print " the total number of bases in DNA is $DNA_num_bases\n";
# percentage of GC
my $GC=($contig=~s/GC/GC/g); 
my $GC_percent={($GC/$Tot_DNA}*100;
print "the percentage of GC is $GC_percent\n";
# percentage of AT
my $AT_percent={($T+$A)/$Tot_DNA)*100;
print "The percentage of AT is $AT_percent\n";
