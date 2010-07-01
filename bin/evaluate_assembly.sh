#!/bin/bash
#
# Daniel Nilsson, 2010
#
# daniel.nilsson@ki.se, daniel.k.nilsson@gmail.com
#
# USAGE: evaluate_assembly.sh assembly_contifs_fasta
#
# Expects: * reference sequences for bmWb and 
#          * blast+ binaries installed in $BLASTBINDIR
#          * several of the other dimma-scripts in $BINDIR
#
# Environment: Several bash environment variables influence the pipeline behaviour
#            BINDIR           (./bin)
#            BLASTBINDIR      (~/src/ncbi-blast-2.2.23+/bin)
#            EXONERATEBINDIR  (~/src/exonerate-2.2.0-x86_64/bin)
#            forceupdate      (no)
#
# Design note: output to log should not be under time-stamp control
#
# check input environment variables and set unset ones to default values

if [ -z "$BLASTBINDIR" ]
then
    BLASTBINDIR=~/src/ncbi-blast-2.2.23+/bin
fi

if [ -z "$BINDIR" ]
then
    BINDIR=~/di/bin
fi

if [ -z "$EXONERATEBINDIR" ]
then   	
    EXONERATEBINDIR=~/src/exonerate-2.2.0-x86_64/bin
fi

if [ -z "$AUGUSTUSBIN" ]
then   	
    AUGUSTUSBIN=~/src/augustus/bin/augustus
fi

if [ -z "$HMMERBINDIR" ]
then
    HMMERBINDIR=~/src/hmmer-3.0/src
fi

if [ -z "$PHRAPBIN" ]
then
    PHRAPBIN=~/src/phrap/phrap
fi

if [ -z "$DBDIR" ]
then
    DBDIR=~/db
fi

# set forceupdate=yes to run all available analyses, even if the file modification times advise against it
if [ -z "$forceupdate" ]
then
    forceupdate=no
fi

# set desired number of concurrent processes
# prefer an already set value of $NPROC, or use nproc to return it if available
if [ -z "$NPROC" ]
then
    NPROCBIN=`which nproc`
    if [ -x $NPROCBIN ] 
    then
	NPROC=`$NPROCBIN`
    fi

    if [ -z "$NPROC" ] 
    then 
	NPROC=1
    fi
fi

function needsUpdate()
{
    # USAGE: needsUpdate(target, prereq [, prereq]*)
    # return true (needsupdate=yes) if target does not yet exist, is older than its prereqs or forceupdate=yes is in effect.

    needsupdate="no"
    
    if [ "$forceupdate" = "yes" ] 
    then
	needsupdate="yes"
    fi

    target=$1;
    
    for prereq in ${@:2}
    do
	if [ $target -ot $prereq ]
	then
	    needsupdate="yes"
	fi
    done
    
    [ "$needsupdate" = "yes" ]
}

# command line parameters 
if [ $# -lt 1 ]
then
        echo "USAGE: ${0##*/} assembly_contigs_fasta_file"
        exit 1
fi

contigs=$1
log=$contigs.eval.log
cat /dev/null > $log

# remove contaminants 

# add actual lengths to name, legacy for the % of contig length criteria evaluation
contigsclen=${contigs}.clen
if needsUpdate $contigsclen $contigs
then
    $BINDIR/add_length_to_fasta_name.pl < $contigs > ${contigs}.clen
fi

# generalise to any input species 
# though other species may need other Evals a.s.o.

function separatespecies()
{
    mycontigs=$1
    ref=$2
    species=$3

    if needsUpdate ${ref}.nsq $ref $BLASTBINDIR/makeblastdb
    then
	$BLASTBINDIR/makeblastdb -in $ref -dbtype nucl
    fi
    
    mycontigswithout=${mycontigs}.no${species}

    mycontigswith=${mycontigs}.${species}

    mycontigsbln=${mycontigs}.${species}.of6.bln

    if needsUpdate $mycontigsbln $mycontigs $BLASTBINDIR/blastn
    then
	$BLASTBINDIR/blastn -num_threads $NPROC -db $ref -query $mycontigs -outfmt 6 -out $mycontigsbln
    fi

    mycontigscf=$mycontigs.${species}.contaminant_filter
    if needsUpdate $mycontigscf $mycontigsbln $BINDIR/catch_contaminants.pl
    then
	$BINDIR/catch_contaminants.pl < $mycontigsbln > $mycontigscf

	cut -f1 $mycontigscf > ${mycontigscf}.ids
    fi
    
    echo -n "${species} screen found: " >> $log
    cat $mycontigscf |sort -k4,4n |awk 'BEGIN {num = 0 ; sum =0 } {sum = sum +$4 ; num =num+1} END { print num " contigs of total length " sum; }' >> $log

    if needsUpdate $mycontigswithout $mycontigscf $EXONERATEBINDIR/fastaremove
    then
	$EXONERATEBINDIR/fastaremove $mycontigs ${mycontigscf}.ids > $mycontigswithout    
    fi

    if needsUpdate $mycontigswith $mycontigscf $EXONERATEBINDIR/fastafetch
    then
	if needsUpdate ${mycontigs}.index $mycontigs $EXONERATEBINDIR/fastaindex
	then
	    $EXONERATEBINDIR/fastaindex -f $mycontigs -i ${mycontigs}.index
	fi

	$EXONERATEBINDIR/fastafetch -f $mycontigs -F -q ${mycontigscf}.ids -i ${mycontigs}.index > $mycontigswith
    fi        
}

# split endosymbiont contigs from nuclear genome 

mycontigs=${contigsclen}
ref=NC_006833.fasta
species=dWb

separatespecies $mycontigs $ref $species

# split mitochondrial contigs from nuclear genome 
mycontigs=${contigsclen}.nodWb
ref=NC_005305_dimmi_mito.fasta
species=mito

separatespecies $mycontigs $ref $species
nuclearcontigs=${contigsclen}.nodWb.nomito

# evaluate EST coverage
covertest=dimmitis_ests_entrez_2009-09-29.fasta

# cluster..
clusters=${covertest}.clusters
if needsUpdate $clusters $covertest $PHRAPBIN
then
    $PHRAPBIN $covertest -ace

    cat ${covertest}.contigs ${covertest}.singlets > $clusters
fi

covertestclen=${clusters}.clen
if needsUpdate $covertestclen $clusters
then
    $BINDIR/add_length_to_fasta_name.pl < $clusters > $covertestclen
fi

if needsUpdate ${nuclearcontigs}.nsq $ref $BLASTBINDIR/makeblastdb
then
    $BLASTBINDIR/makeblastdb -in $nuclearcontigs -dbtype nucl
fi

guidecover=$covertestclen.$nuclearcontigs.guidecover
covertestclencontigsbln=${covertestclen}.${nuclearcontigs}.of6.bln

if needsUpdate $covertestclencontigsbln $mycontigs $BLASTBINDIR/blastn
then
    $BLASTBINDIR/blastn -num_threads $NPROC -db $nuclearcontigs -query $covertestclen -outfmt 6 -out $covertestclencontigsbln
fi

if needsUpdate $guidecover $covertestclencontigsbln $BINDIR/check_guide_coverage.pl
then
    $BINDIR/check_guide_coverage.pl < $covertestclencontigsbln > $guidecover
fi

total=`grep -c \> $covertest`
totalclust=`grep -c \> $clusters`
echo "Results from cover test with $covertest ($total sequences, $totalclust clusters)." >> $log

foundok=`grep -c FOUND_OK $guidecover`
fragmented=`grep FOUND_OK $guidecover |awk 'BEGIN { num=0} ($5>1) {num=num+1} END {print num}'`
unfragmented=`grep FOUND_OK $guidecover |awk 'BEGIN { num=0} ($5==1) {num=num+1} END {print num}'`
echo "Found ok: "$foundok" United: "$unfragmented" Fragmented: "$fragmented >> $log

hits=`grep -c HITS_ONLY $guidecover`
hitsfragmented=`grep HITS_ONLY $guidecover |awk 'BEGIN { num=0} ($5>1) {num=num+1} END {print num}'`
hitsunfragmented=`grep HITS_ONLY $guidecover |awk 'BEGIN { num=0} ($5==1) {num=num+1} END {print num}'`
echo "Found hits: "$hits" United: "$hitsunfragmented" Fragmented: "$hitsfragmented >> $log

nothit=`grep -c NOT_FOUND $guidecover`
echo "Not accepted hits: $nothit" >> $log

echo "Total not found: "$(( $totalclust - $foundok )) >> $log

# Evaluate gene structure and contents

# Augustus gene prediction

# gene finding

model=brugia
augustusgff=${nuclearcontigs}.augustus.$model.gff3

if needsUpdate $augustusgff $nuclearcontigs $AUGUSTUSBIN 
then
    echo -n "Start Augustus..."
    $AUGUSTUSBIN --species=$model --outfile=$augustusgff $nuclearcontigs --gff3=on

    echo "Augustus Done."
fi

grep -v \# $augustusgff |awk 'BEGIN {num=0} ($3=="gene") { num=num+1 } END {print "Augugstus found " num " genes"}' >> $log

augustuspep=${augustusgff%%gff3}pep
if needsUpdate $augustuspep $augustusgff $BINDIR/get_protein_from_augustus_gff.pl
then
    $BINDIR/get_protein_from_augustus_gff.pl < $augustusgff > $augustuspep
fi

augustusembl=${augustusgff%%gff3}embl
if needsUpdate $augustusembl $augustusgff $nuclearcontigs $BINDIR/merge_fasta_and_gff.pl
then
    $BINDIR/merge_fasta_and_gff.pl -g $augustusgff -f $nuclearcontigs -o $augustusembl -t embl
fi

genestructurestats=$augustusembl.stats
if needsUpdate $genestructurestats $augustusembl $BINDIR/gene_structure_stats.pl 
then
    $BINDIR/gene_structure_stats.pl $augustusembl > $genestructurestats
fi

cat $genestructurestats >> $log

# GLIMMER on the wolbachia..

# Then...
# hmmer search,
# this fits better in a real annotation pipeline.. but, I guess this will become the annotation pipeline. 
# in that case, running LGIC finding/tree construction is rather interesting as well! this also needs hmmer3-filter

hmmerout=$augustuspep.pfamA.tc.hmmer_out
if needsUpdate $hmmerout $augustuspep $HMMERBINDIR/hmmsearch
then
    $HMMERBINDIR/hmmsearch --cut_tc -o $hmmerout --cpu $NPROC $DBDIR/Pfam-A.hmm $augustuspep
fi

# extract fasta files with an id
#*manual cut & paste rectangles from hmmer-PFAM results*
#cat 2009-12-07_GEU-1_velvet-h31_SCF.augustus_brugia.Neur_chan_LBD.ids 2009-12-07_GEU-1_velvet-h31_SCF.augustus_brugia.Neur_chan_memb.ids |sort |uniq|sort >2009-12-07_GEU-1_velvet-h31_SCF.augustus_brugia.Neur_chan_LBD_and_memb.ids
#../src/exonerate-2.2.0-x86_64/bin/fastaindex -f 2009-12-07_GEU-1_velvet-h31_SCF.augustus_brugia.pep -i 2009-12-07_GEU-1_velvet-h31_SCF.augustus_brugia.index
#../src/exonerate-2.2.0-x86_64/bin/fastafetch -f 2009-12-07_GEU-1_velvet-h31_SCF.augustus_brugia.pep -F -q 2009-12-07_GEU-1_velvet-h31_SCF.augustus_brugia.Neur_chan_LBD_and_memb.ids -i 2009-12-07_GEU-1_velvet-h31_SCF.augustus_brugia.index > 2009-12-07_GEU-1_velvet-h31_SCF.augustus_brugia.Neur_chan_LBD_and_memb.pep

