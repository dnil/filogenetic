#!/bin/bash
#
# (c) Daniel Nilsson, 2010
#
# daniel.nilsson@ki.se, daniel.k.nilsson@gmail.com
#
# Package released under the Perl Artistic License.
#
# USAGE: evaluate_assembly.sh assembly_contigs_fasta
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

: <<'POD_INIT'

=head1 NAME

evaluate_assembly.sh - evaluate D. immitis assembly completeness

=head1 AUTHOR

Daniel Nilsson, daniel.nilsson@izb.unibe.ch, daniel.nilsson@ki.se, daniel.k.nilsson@gmail.com

=head1 LICENSE AND COPYRIGHT

Copyright 2009, 2010 held by Daniel Nilsson. The package is released for use under the Perl Artistic License.

=head1 SYNOPSIS

USAGE: C<evaluate_assembly.sh contigs.fa>

=head1 DESCRIPTION

Evaluates a filarial worm assembly based on EST completeness coverage stats.
Also checks for Wolbachia endosymbiont and mitochondrial DNA, calls tentatively calls 
genes using Augustus, runs a Hmmer3 PFAM screen on the conceptual peptides and checks motif content,
especially with regard to LGICs.

Uses pipelinefunk.sh - see this for more usage info.

=head1 DEPENDENCIES

Expects:

=over 4

=item * 

Assembly fasta file.

=item *
          
blast+, exonerate, Augustus, Hmmer3 and PHRAP installed. Separate ENV environment variables should be set to point to the installations of these. See below.

=item *

BioPerl for gene structure stats.

=item *

PFAM, version suitable for Hmmer3 (e.g. v24).

=item *

Wolbachia genome, DI mito genome, DI ESTs - a version of these should be available in the data dir.

=back

=head1 SHELL ENVIRONMENT VARIABLES

Several environment variables influence the pipeline behaviour.

=over 4

=item BINDIR [path (~/sandbox/filogenetic/bin)]          

Directory where the rest of the tools live.

=item BLASTBINDIR [path (~/src/ncbi-blast-2.2.23+/bin)] 

Directory where the NCBI blast+ binaries reside.

=item EXONERATEBINDIR [path (~/src/exonerate-2.2.0-x86_64/bin)]

Directory where the Exonerate binaries reside.

=item HMMERBINDIR [path (~/src/hmmer-3.0/bin)]

Directory where the Hmmer-3 binaries reside.

=item AUGUSTUSBIN [path (~/src/augustus/bin/augustus)]

Path to the Augustus binary.

=item PHRAPBIN [path (~/src/phrap/phrap)]

Path to the Augustus binary.

=item DBDIR [path (~/db)]

Path to where the databases used can be found. Currently only uses PFAM.

=back

A few more are inherited from pipelinefunk.sh: DIRECTIVE, forceupdate

POD_INIT

if [ -z "$BLASTBINDIR" ]
then
    BLASTBINDIR=~/src/ncbi-blast-2.2.23+/bin
fi

if [ -z "$BINDIR" ]
then
    BINDIR=~/sandbox/filogenetic/bin
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

# For the use of any local perl modules (none for this package yet!)
#export PERL5LIB=$PERL5LIB:$BINDIR

# uses pipelinefunk.sh for needsUpdate, registerFile, NPROC etc.

. $BINDIR/pipelinefunk.sh

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
registerFile $contigsclen temp
if needsUpdate $contigsclen $contigs
then
    $BINDIR/add_length_to_fasta_name.pl < $contigs > ${contigs}.clen
fi

# generalise to any input species 
# though other species may need other Evals a.s.o.
: << 'FUNCTION_DOC'

=head1 FUNCTIONS

=head2 separatespecies(mycontigs, ref, species)

Pull apart separate species from the assembly, e.g. mito DNA or endosymbiont DNA.

FUNCTION_DOC

function separatespecies()
{
    mycontigs=$1
    ref=$2
    species=$3

    if [ ! -e $ref ]
    then
	echo Could not find $ref. Testing other locations.. 1>&2
	echo How about ${0%%bin*}data/$ref or ${0%%bin*}data/$ref.bz2? 1>&2
	if [ -e ${0%%bin*}data/$ref ]
	then
	    ln -s ${0%%bin*}data/$ref .
	elif [ -e ${0%%bin*}data/$ref.bz2 ] 
	then
	    bunzip2 ${0%%bin*}data/$ref.bz2
	    ln -s ${0%%bin*}data/$ref .
	else
	    echo "OOPS! Reference file $ref not found!"
	fi
    fi

    registerFile $ref.nsq temp
    registerFile $ref.nin temp
    registerFile $ref.nhr temp

    if needsUpdate ${ref}.nsq $ref $BLASTBINDIR/makeblastdb
    then
	$BLASTBINDIR/makeblastdb -in $ref -dbtype nucl
    fi

    mycontigswithout=${mycontigs}.no${species}

    mycontigswith=${mycontigs}.${species}

    mycontigsbln=${mycontigs}.${species}.of6.bln

    registerFile $mycontigsbln temp
    if needsUpdate $mycontigsbln $mycontigs $BLASTBINDIR/blastn
    then
	$BLASTBINDIR/blastn -num_threads $NPROC -db $ref -query $mycontigs -outfmt 6 -out $mycontigsbln
    fi

    mycontigscf=$mycontigs.${species}.contaminant_filter

    registerFile $mycontigscf temp
    registerFile $mycontigscf.ids temp
    if needsUpdate $mycontigscf $mycontigsbln $BINDIR/catch_contaminants.pl
    then
	$BINDIR/catch_contaminants.pl < $mycontigsbln > $mycontigscf

	cut -f1 $mycontigscf > ${mycontigscf}.ids
    fi
    
    echo -n "${species} screen found: " >> $log
    cat $mycontigscf |sort -k4,4n |awk 'BEGIN {num = 0 ; sum =0 } {sum = sum +$4 ; num =num+1} END { print num " contigs of total length " sum; }' >> $log

    registerFile $mycontigswithout result
    if needsUpdate $mycontigswithout $mycontigscf $EXONERATEBINDIR/fastaremove
    then
	$EXONERATEBINDIR/fastaremove $mycontigs ${mycontigscf}.ids > $mycontigswithout    
    fi

    registerFile $mycontigswith result
    if needsUpdate $mycontigswith $mycontigscf $EXONERATEBINDIR/fastafetch
    then
	if needsUpdate ${mycontigs}.index $mycontigs $EXONERATEBINDIR/fastaindex
	then
	    if [ -e $mycontigs.index ]
	    then
		rm ${mycontigs}.index
	    fi
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
if [ ! -e $covertest ] 
then
    echo Could not find $covertest. Testing other locations..  >&2 
    # try if it is available in data?
    if [ -e ${0%%bin*}data/$covertest ] 
    then
	ln -s ${0%%bin*}data/$covertest .
    elif [ -e ${0%%bin*}data/$covertest.bz2 ]
    then
	bunzip2 ${0%%bin*}data/$covertest.bz2
	ln -s ${0%%bin*}data/$covertest .
    else
	echo "OOPS! Covertest file $covertest not found!"
    fi
fi

# cluster..
clusters=${covertest}.clusters
registerFile $clusters temp
if needsUpdate $clusters $covertest $PHRAPBIN
then
    $PHRAPBIN $covertest -ace

    cat ${covertest}.contigs ${covertest}.singlets > $clusters
fi

covertestclen=${clusters}.clen
registerFile $covertestclen temp
if needsUpdate $covertestclen $clusters
then
    $BINDIR/add_length_to_fasta_name.pl < $clusters > $covertestclen
fi

registerFile ${nuclearcontigs}.nsq temp
registerFile ${nuclearcontigs}.nin temp
registerFile ${nuclearcontigs}.nhr temp
if needsUpdate ${nuclearcontigs}.nsq $ref $BLASTBINDIR/makeblastdb
then
    $BLASTBINDIR/makeblastdb -in $nuclearcontigs -dbtype nucl
fi

guidecover=$covertestclen.$nuclearcontigs.guidecover
covertestclencontigsbln=${covertestclen}.${nuclearcontigs}.of6.bln

registerFile $covertestclencontigsbln temp 
if needsUpdate $covertestclencontigsbln $mycontigs $BLASTBINDIR/blastn
then
    $BLASTBINDIR/blastn -num_threads $NPROC -db $nuclearcontigs -query $covertestclen -outfmt 6 -out $covertestclencontigsbln
fi

registerFile $guidecover temp
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

registerFile $augustusgff result
if needsUpdate $augustusgff $nuclearcontigs $AUGUSTUSBIN 
then
    echo -n "Start Augustus..."
    $AUGUSTUSBIN --species=$model --outfile=$augustusgff $nuclearcontigs --gff3=on

    echo "Augustus Done."
fi

grep -v \# $augustusgff |awk 'BEGIN {num=0} ($3=="gene") { num=num+1 } END {print "Augugstus found " num " genes"}' >> $log


augustuspep=${augustusgff%%gff3}pep
registerFile $augustuspep result
if needsUpdate $augustuspep $augustusgff $BINDIR/get_protein_from_augustus_gff.pl
then
    $BINDIR/get_protein_from_augustus_gff.pl < $augustusgff > $augustuspep
fi

augustusembl=${augustusgff%%gff3}embl
registerFile $augustusembl temp
if needsUpdate $augustusembl $augustusgff $nuclearcontigs $BINDIR/merge_fasta_and_gff.pl
then
    $BINDIR/merge_fasta_and_gff.pl -g $augustusgff -f $nuclearcontigs -o $augustusembl -t embl
fi

genestructurestats=$augustusembl.stats
registerFile $genestructurestats temp
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
registerFile $hmmerout result
if needsUpdate $hmmerout $augustuspep $HMMERBINDIR/hmmsearch
then
    $HMMERBINDIR/hmmsearch --cut_tc -o $hmmerout --cpu $NPROC $DBDIR/Pfam-A.hmm $augustuspep
fi

hmmer_all_ids=${hmmerout%%.hmmer_out}.allids
hmmer_neur_chan_ids=${hmmerout%%.hmmer_out}.neurchan.ids
hmmer_neur_chan_pep=${hmmerout%%.hmmer_out}.neurchan.pep

registerFile $hmmer_all_ids temp
registerFile $hmmer_neur_chan_ids temp
registerFile $hmmer_neur_chan_pep result

if needsUpdate $hmmer_neur_chan_pep $hmmerout $BINDIR/get_hmmer3search_domain_coords.pl
then
    $BINDIR/get_hmmer3search_domain_coords.pl $hmmerout > $hmmer_all_ids
    grep Neur_chan $hmmer_all_ids|cut -f 2  > $hmmer_neur_chan_ids.tmp
    grep Lig_chan $hmmer_all_ids|cut -f 2  >> $hmmer_neur_chan_ids.tmp
    sort $hmmer_neur_chan_ids.tmp | uniq > $hmmer_neur_chan_ids
    
    if [ -e ${augustuspep}.index ] 
    then
	rm ${augustuspep}.index
    fi

    $EXONERATEBINDIR/fastaindex -f $augustuspep -i ${augustuspep}.index 
    $EXONERATEBINDIR/fastafetch -f $augustuspep -F -q $hmmer_neur_chan_ids -i ${augustuspep}.index > ${hmmer_neur_chan_pep}

fi
