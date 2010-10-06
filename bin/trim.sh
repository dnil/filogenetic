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

USAGE: C<trim.sh mate1.fastq mate2.fastq>

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

POD_INIT

if [ -z "$FASTXBINDIR" ]
then
    BLASTBINDIR=~/src/fastx-toolkit/bin
fi

if [ -z "$BINDIR" ]
then
    BINDIR=~/sandbox/filogenetic/bin
fi

if [ -z "$VELVETBINDIR" ]
then   	
    VELVETBINDIR=~/src/velvet
fi

quallim=20

# For the use of any local perl modules (none for this package yet!)
#export PERL5LIB=$PERL5LIB:$BINDIR

# uses pipelinefunk.sh for needsUpdate, registerFile, NPROC etc.

. $BINDIR/pipelinefunk.sh

# command line parameters 
if [ $# -lt 2 ]
then
        echo "USAGE: ${0##*/} mate1.fastq mate2.fastq"
        exit 1
fi

mate1=$1
mate2=$2
log=$mate1.eval.log
cat /dev/null > $log


#shuffle 
velvetshuffle=$mate1.velvetshuffle.fastq
registerFile $velvetshuffle result
if needsUpdate $velvetshuffle $mate1 $mate2 $VELVETBINDIR/shuffleSequences_fastq.pl
then
    $VELVETBINDIR/shuffleSequences_fastq.pl $mate1 $mate2 $velvetshuffle
fi

#endtrim
endtrim=${velvetshuffle%%.fastq}.endtrimq${quallim}.fastq
registerFile $entrim temp
if needsUpdate $endtrim $velvetshuffle
then
    $FASTXBINDIR/fastq_quality_trimmer -t $quallim -l 0 -i $velvetshuffle -o $endtrim
fi

reapN=${endtrim%%.fastq}.noN.fastq
registerFile $reapN temp
if needsUpdate $reapN $endtrim $BINDIR/qualreap.pl
then
    $BINDIR/qualreap.pl < $endtrim > $reapN
fi

pairsonly=${reapN%%.fastq}.pairsonly.fastq
registerFile $pairsonly temp
if needsUpdate $pairsonly $reapN $BINDIR/keeponlypairs.pl
then
    $BINDIR/keeponlypairs.pl < $reapN > $pairsonly
fi

longenough=${pairsonly%%.fastq}.longenough.fastq
registerFile $longenough result
if needsUpdate $longenough $pairsonly $BINDIR/killpairifshort.pl
then
    $BINDIR/killpairifshort.pl < $pairsonly > $longenough
fi

# ../src/fastx-toolkit/bin/fastx_reverse_complement -i GEU-7b.velvetshuffle.endtrimq20.noN.fastq -o GEU-7b.velvetshuffle.endtrimq20.noN.revcomp.fastq

# STDERR to log?

