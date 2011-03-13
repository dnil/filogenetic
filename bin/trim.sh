#!/bin/bash
#
# (c) Daniel Nilsson, 2010
#
# daniel.nilsson@ki.se, daniel.k.nilsson@gmail.com
#
# Package released under the Perl Artistic License.
#
# USAGE: trim.sh mate1.fastq mate2.fastq
#
# Expects: * fastx-tools 
#          * velvet (velvet_shuffle.pl)
#          * several of the other package scripts in $BINDIR
#
# Environment: Several bash environment variables influence the pipeline behaviour
#            BINDIR           (./bin)
#            FASTXBINDIR  (~/src/fastx-toolkit/bin)
#            forceupdate      (no)
#
# Design note: output to log should not be under time-stamp control
#
# check input environment variables and set unset ones to default values

: <<'POD_INIT'

=head1 NAME

trim.sh - trim reads to use in D. immitis assembly 

=head1 AUTHOR

Daniel Nilsson, daniel.nilsson@scilifelab.se, daniel.nilsson@ki.se, daniel.k.nilsson@gmail.com

=head1 LICENSE AND COPYRIGHT

Copyright 2009, 2010 held by Daniel Nilsson. The package is released for use under the Perl Artistic License.

=head1 SYNOPSIS

USAGE: C<trim.sh mate1.fastq mate2.fastq [outbase]>

=head1 DESCRIPTION

Trim reads.

Uses pipelinefunk.sh - see this for more usage info.

=head1 DEPENDENCIES

Expects:

=over 4

=item * 

Reads in mate separated files.

=item * 

fastx-toolkit

=item *

velvet (for velvet_shuffle.pl)

=back

=head1 SHELL ENVIRONMENT VARIABLES

Several environment variables influence the pipeline behaviour.

=over 4

=item BINDIR [path (~/sandbox/filogenetic/bin)]          

Directory where the rest of the tools live.

=item FASTXBINDIR [path (~/src/fastx-toolkit/bin)]

Directory where the Exonerate binaries reside.

=item VELVETBINDIR [path (~/src/velvet]

Directory where the velvet binaries reside.

=back

A few more are inherited from pipelinefunk.sh: DIRECTIVE, forceupdate.

=cut

POD_INIT

if [ -z "$FASTXBINDIR" ]
then
    FASTXBINDIR=~/src/fastx-toolkit/bin
fi

if [ -z "$BINDIR" ]
then
    BINDIR=~/sandbox/filogenetic/bin
fi

if [ -z "$VELVETBINDIR" ]
then   	
    VELVETBINDIR=~/src/velvet
fi

if [ -z "$truncate" ]
then 
    truncate="no"
fi

if [ -z "$trunclen" ]
then
    trunclen=50
fi

if [ -z "$revcomp" ]
then
    revcomp="no"
fi

if [ -z "$quallim" ]
then
    quallim=20
fi

# For the use of any local perl modules (none for this package yet!)
#export PERL5LIB=$PERL5LIB:$BINDIR

# uses pipelinefunk.sh for needsUpdate, registerFile, NPROC etc.

. $BINDIR/pipelinefunk.sh

# command line parameters 
if [ $# -lt 2 ]
then
        echo "USAGE: ${0##*/} mate1.fastq mate2.fastq [newoutbasename]"
        exit 1
fi

mate1=$1
mate2=$2

outname=$mate1
if [ $# -eq 3 ]
then
    outname=$3
fi

log=$outname.log
#cat /dev/null > $log
rundate=`date`
echo "$0 $mate1 $mate2 $outname revcomp=$revcomp truncate=$truncate trunclen=$trunclen quallim=$quallim pwd=$PWD ($rundate)" >> $log

updates="no"

#shuffle 
velvetshuffle=$outname.velvetshuffle.fastq
registerFile $velvetshuffle result
if needsUpdate $velvetshuffle $mate1 $mate2 $VELVETBINDIR/shuffleSequences_fastq.pl
then
    $VELVETBINDIR/shuffleSequences_fastq.pl $mate1 $mate2 $velvetshuffle 2>> $log
    updates="yes"
fi

#truncate?
if [ "$truncate" = "yes" ]
then
    truncated=${velvetshuffle%%.fastq}.trunc${trunclen}.fastq
    registerFile $truncated temp
    if needsUpdate $truncated $velvetshuffle
    then
	echo "$FASTXBINDIR/fastx_trimmer -l $trunclen -i $velvetshuffle -o $truncated" >> $log
	$FASTXBINDIR/fastx_trimmer -l $trunclen -i $velvetshuffle -o $truncated 2>> $log
	updates="yes"
    fi
else
    truncated=$velevetshuffle
fi

#endtrim
endtrim=${truncated%%.fastq}.endtrimq${quallim}.fastq
registerFile $endtrim temp
if needsUpdate $endtrim $truncated
then
    echo "$FASTXBINDIR/fastq_quality_trimmer -t $quallim -l 0 -i $truncated -o $endtrim" >> $log
    $FASTXBINDIR/fastq_quality_trimmer -t $quallim -l 0 -i $truncated -o $endtrim 2>> $log
    updates="yes"
fi

# remove N containing reads
reapN=${endtrim%%.fastq}.noN.fastq
registerFile $reapN temp
if needsUpdate $reapN $endtrim $BINDIR/qualreap.pl
then
    echo "$BINDIR/qualreap.pl < $endtrim > $reapN" >> $log
    $BINDIR/qualreap.pl < $endtrim > $reapN 2>> $log
    updates="yes"
fi

# retain only even pairs
pairsonly=${reapN%%.fastq}.pairsonly.fastq
registerFile $pairsonly temp
if needsUpdate $pairsonly $reapN $BINDIR/keeponlypairs.pl
then
    echo "$BINDIR/keeponlypairs.pl < $reapN > $pairsonly" >>$log
    $BINDIR/keeponlypairs.pl < $reapN > $pairsonly 2>> $log
    updates="yes"
fi

# get rid of pairs with a too short mate
longenough=${pairsonly%%.fastq}.longenough.fastq
registerFile $longenough result
if needsUpdate $longenough $pairsonly $BINDIR/killpairifshort.pl
then
    echo "$BINDIR/killpairifshort.pl < $pairsonly > $longenough" >> $log
    $BINDIR/killpairifshort.pl < $pairsonly > $longenough 2>> $log
fi

if [ "$revcomp" = "yes" ] 
then 
    revcomped=${longenough%%.fastq}.revcomp.fastq
    registerFile $revcomp result

    if needsUpdate $revcomped $longenough
    then
	echo "$FASTXBINDIR/fastx_reverse_complement -i $longenough -o $revcomped" >>$log
	$FASTXBINDIR/fastx_reverse_complement -i $longenough -o $revcomped 2>> $log
	updates="yes"
    fi
fi

#
if [ "$updates" = "yes" ]
then
    final=$longenough 
    if [ "$revcomp" = "yes" ]
    then
	final=$revcomped
    fi

    echo "Number of incorrect pairs (should always be 0!): " >> $log
    grep "^@" $final |cut -d\# -f 1 |uniq -c |grep -cv \ 2\  >> $log
    echo "Total number of reads after trimming: " >> $log
    grep -c "^@" $final >> $log
    echo "Total number of bases left after trimming: " >> $log
    
    perl -e 'my $total=0; my $n=0; while (<STDIN>) { if ($n ==1) { chomp; $total+=length ; $n++;} elsif ($n == 3) {$n=0 } else {$n++}} print "Found $total bp.\n";' < $final >> $log

    completion=`date`
    echo "($completion) Project brought up to date." >> $log
else 
    echo "Project was already up to date." >> $log
fi

# STDERR to log?

echo "Done. See log $log for details."
