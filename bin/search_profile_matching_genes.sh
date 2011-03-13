#!/bin/bash

# usage:
# search_profile_matching_genes.sh [PFAM_motif] [NOTREES]
#

#
# DOC/TODO: 
# * script wrap db split and blast db construction
# * stand alone motif extraction (from pep seq, that is) (package level todo..)
#
# * introduce pipelinefunk wrapper, and deconvolute some obscure passages

:<<DOC

=head1 NAME

search_profile_matching_genes.sh - find genes with defined peptide motif in draft nucleotide sequence

=head1 AUTHOR

Daniel Nilsson, daniel.nilsson@izb.unibe.ch, daniel.nilsson@ki.se, daniel.k.nilsson@gmail.com

=head1 LICENSE AND COPYRIGHT

Copyright 2009-2011 held by Daniel Nilsson. The package is realesed for use under the Perl Artistic License.

=head1 SYNOPSIS

USAGE: C<search_profile_matching_genes.sh [PFAM_motif] [NOTREES]>

=head1 DESCRIPTION

Find sequenes bearing a given motif in genomic draft seqeucence. Uses
GeneWise after some non-stringent blast filtering of contigs (split
into performance-wise convenient pieces), well knowing the author of
GeneWise might not have approved. Classifies the predicted peptides by
means of phylogenetic tree construction over the extracted motif with
a reference seed set.

=head1 DEPENDENCIES

NCBI blast, hmmsearch (v3), genewise, clustalw and EMBOSS on the search PATH.

blast+ is on the todolist, not really for speed but
avoiding legacy versions would make installation on modern systems
more straightforward). Exonerate would be an optional replacement for
GeneWise/Wise2, but at the time of implementation the gene modelling
part was not really that advanced yet, and there was no support for
pHMM input; so it would be seed protein only.

=head1 ENVIRONMENT VARIABLES

A number of environment variables affect the run. Note in particular the 
$SEEDQUERYFILE, which is important for classification.

=over 4

=item MYHOME [directory ($HOME/nematodes)]

Main project directory, typically containing a number of subdirectories.

=item DBDIR [directory ($MYHOME/db)]

Database directory, containing fasta files of genomes of interest,
truncated into e.g. 100kb fragments to keep processing time down in
the GeneWise search.

=item SEQTEMP [directory ($MYHOME/seq_temp)]

Tempdir used for sequences etc.

=item SEARCHES [directory ($MYHOME/searches)]

Results (and some intermediate files).

=item SEEDQUERYFILE [peptidefastafile (seed_ce_go0005230_uniprothumanlgic_LRseqHcDEG3s.pep)]

Seed sequences used in blast screen and for orientation in the final
alignment/tree construction.

=item DBS [fasta file list (cut -f1 $DBDIR/db_initials)]

List of the fasta databases to use.

N.B. that the $DBDIR/db_initials file needs to be populated anyway to
have nice short names on the tree!

=item MOTIF [pHMM-file (Neur_chan_LBD_ls.hmm)]

Profile HMM file readable by GeneWise, e.g. from PFAM, or your own
home brew generated with HMMER.

=item BLASTBINDIR [path (~/src/blast-2.2.24/bin)]

NCBI BLAST legacy version, binary dir.

=item HMMERBINDIR [path (~/src/hmmer-3.0/src)]

Path to Hmmer3 binary dir.

=back

=cut

DOC

# environment variables

if [ -z "$MYHOME" ]
then
    MYHOME=$HOME/nematodes
fi

if [ -z "$DBDIR" ]
then
    DBDIR=$MYHOME/db
fi

if [ -z "$BINDIR" ]
then
    BINDIR=$MYHOME/bin
fi

if [ -z "$SEQTEMP" ]
then
    SEQTEMP=$MYHOME/seq_temp
fi

if [ -z "$SEARCHES" ]
then
    SEARCHES=$MYHOME/searches
fi

if [ -z "$SEEDQUERYFILE" ]
then
    SEEDQUERYFILE=seed_ce_go0005230_uniprothumanlgic_LRseqHcDEG3s.pep
fi

#DBS=cjaponica_4.0.1_supercontigs.100k.L2k.split.fasta

#DBS="c_elegans.WS194.dna.100kOL2k.clustalnames"

if [ -z "$DBS" ]
then 
#    PLATYHELMDBS="schmidtea_mediterranea_31_contigs.100k.L2k.split.fasta sma_v3.1.fasta E_multilocularis_contigs_230108.fa"
#    NEMDBS="Heterorhabditis_bacteriophora-1.2.1.contigs.fa cremanei_contigs.100k.L2k.split.fasta combined_worms_supercontigs200808.fasta b_malayi.WS194.100k.L2k.split.fasta meloidogyne_entrez.fasta Mhapla_wgs080929.wgs.100k.L2k.split.fasta pristionchus_pacificus_contigs.fa t_spiralis1.0-contigs.fa cjaponica_4.0.1_supercontigs.100k.L2k.split.fasta Cbrenneri_PB2801_6.0.1.100k.L2k.split.fasta c_briggsae.100kOL2k.fasta S_ratti_contigs.141008 Ascaris_suum.WGS.contigs nembase3.clustalnames nematode_net_4_species_est_clusters.fasta N_brasiliensis_454_contigs_033109.fasta TCIR.supercontigs.Phusion.120209.fasta"
#    VERTDBS="tfru.100kL2k.fasta btu_mdass_full.100k.L2k.split.fasta Homo_sapiens.NCBI36.54.dna.toplevel.100k.L2k.split.fasta Xenopus_tropicalis.JGI4.1.54.dna.100k.L2k.split.fasta Gasterosteus_aculeatus.BROADS1.54.dna.100k.L2k.split.fasta Gallus_gallus.WASHUC2.54.dna.100k.L2k.split.fasta Felis_catus.CAT.54.dna.100k.L2k.split.fasta Danio_rerio.Zv8.54.dna.toplevel.100k.L2k.split.fasta Canis_familiaris.BROADD2.54.dna.100k.L2k.split.fasta"
#    CELDB="c_elegans.WS194.dna.100kOL2k.clustalnames"
#    INSECTDBS="Aedes_aegypti.AaegL1.52.100k.L2k.split.fasta Anopheles_gambiae.AgamP3.52.100k.L2k.split.fasta Culex_quinquefasciatus.CpipJ1.52.100k.L2k.split.fasta Drosophila_melanogaster.BDGP5.4.52.100k.L2k.split.fasta Ixodes_scapularis.IscaW1.52.100k.L2k.split.fasta "
#    plantDBS="Arabidopsis_lyrata.JGI8X.52.100k.L2k.split.fasta Arabidopsis_thaliana.TAIR8.52.100k.L2k.split.fasta Oryza_glaberrima.unspecified.52.100k.L2k.split.fasta Oryza_sativa_indica.Jan_2005.52.100k.L2k.split.fasta Oryza_sativa_japonica.TIGR5.52.100k.L2k.split.fasta Populus_trichocarpa.jgi2004.52.100k.L2k.split.fasta Sorghum_bicolor.JGI2007.52.100k.L2k.split.fasta Vitis_vinifera.8X.52.100k.L2k.split.fasta"

#    DBS=$NEMDBS" "$PLATYHELMDBS" "$VERTDBS" "$INSECTDBS
    DBS=`cut -f1 $DBDIR/db_initials`
fi

# uses pipelinefunk.sh for needsUpdate, registerFile, NPROC etc.

. $BINDIR/pipelinefunk.sh

#
# command line parameters
#

#MOTIF="Neur_chan_LBD_ls.hmm" # Neur_chan_LBD_ls.hmm
#MOTIF="Neur_chan_memb_ls.hmm" # Neur_chan_LBD_ls.hmm
# _ls motifs are "glocal"-search, _fs "swat"-search; prior is most sensitive IF whole motif represented, latter only option for fragments..
# NB: motif file currently required to reside in pwd

if [ "$#" -gt 0 ]
then
    MOTIF=$1
fi

if [ -z "$MOTIF" ] 
then
    MOTIF=Neur_chan_LBD_ls.hmm
fi

if [ -z "$BLASTBINDIR" ]
then 
    BLASTBINDIR=~/src/blast-2.2.24/bin
fi

if [ -z "$HMMERBINDIR" ]
then
    HMMERBINDIR=~/src/hmmer-3.0/src
fi

if [ -z "$WISECONFIGDIR" ]
then
    WISECONFIGDIR=/home/daniel/src/wise2.2.0/wisecfg
fi

DOTREES=1
if [ "$#" -gt 1 ]
then
    STAGE=$2
    if [ "$STAGE" == "NOTREES" ]
    then
	DOTREES=0
    else
	echo Stage description $STAGE not understood. Cowardly bailing out.
	exit 1
    fi
fi

#
# extract motif subsequences from the seed sequence set, for phylogenetically sound alignments
#

seed_hmmsearch=$SEARCHES/${MOTIF}_in_${SEEDQUERYFILE}.hmmsearch
if needsUpdate $seed_hmmsearch $MOTIF $SEEDQUERYFILE
then	
    $HMMERBINDIR/hmmsearch $MOTIF $SEEDQUERYFILE > $seed_hmmsearch
    registerFile $seed_hmmsearch temp
fi

seed_clusteyenames=${SEQTEMP}/${SEEDQUERYFILE}.clusteyenames
if needsUpdate $seed_clusteyenames $SEEDQUERYFILE 
then
    perl -ne 'if (m/^>(\S+)/) { $name = $1; $name=~s/\|/+/g; print ">$name\n"; } else {print;}' < ${SEEDQUERYFILE} > $seed_clusteyenames
    registerFile $seed_clusteyenames temp
fi
    
if needsUpdate todo $seed_hmmsearch $BINDIR/get_hmmer3search_domain_coords.pl
then
    $BINDIR/get_hmmer3search_domain_coords.pl $seed_hmmsearch | awk '{print "seqret '${seed_clusteyenames}':"$1,"'${SEQTEMP}'/"$1"_"$2"-"$3".fasta","-sbegin",$2,"-send",$3;}'|sed -e 's/|/+/g;' > todo
    registerFile todo temp

    chmod a+x todo
    ./todo

    # will create a number of unregistered temp fasta files in seqtemp...
fi

seed_motifs_fasta_clustalnames=${SEQTEMP}/${MOTIF}_in_${SEEDQUERYFILE}.clustalnames
if needsUpdate $seed_motifs_fasta_clustalnames $BINDIR/get_hmmer3search_domain_coords.pl $seed_hmmsearch
then
    seed_hmmsearch_fastalist=${SEQTEMP}/${MOTIF}_in_${SEEDQUERYFILE}.hmmsearch.fastalist	
    $BINDIR/get_hmmer3search_domain_coords.pl $seed_hmmsearch | awk '{print "'$SEQTEMP'/"$1"_"$2"-"$3".fasta"}' |sed -e 's/|/+/g;' > $seed_hmmsearch_fastalist
    registerFile $seed_hmmsearch_fastalist temp

    seed_motifs_fasta=${SEQTEMP}/${MOTIF}_in_${SEEDQUERYFILE}.hmmsearch.fasta
    echo -n > $seed_motifs_fasta
    for file in `cat $seed_hmmsearch_fastalist`; do echo $file | perl -ne '$row = $_; chomp $row; m/(\d+)-(\d+)/; $begin = $1; $end=$2; open FASTA, $row; while (my $fr=<FASTA>) { chomp $fr; if($fr=~m/^>(.+)/) {print $fr."_".$begin."_".$end."\n";} else {print $fr."\n";} } close FASTA;' ; done > $seed_motifs_fasta
	
    registerFile $seed_motifs_fasta temp

    perl -ne 'if (m/^>(\S+)/) { $name = $1; $name=~s/\+/\|/g; print ">$name\n"; } else {print;}' < $seed_motifs_fasta > $seed_motifs_fasta_clustalnames

    registerFile $seed_motifs_fasta_clustalnames result
fi

#
# pre-screen genomes with blast for low stringency tblastn hits to the seed sequences
# then run genewise on the hit contigs only
#

for db in $DBS ; do 
    if needsUpdate $DBDIR/$db.nsq $DBDIR/$db $BLASTBINDIR/formatdb
    then
	$BLASTBINDIR/formatdb -i $DBDIR/$db -p F
	    
	registerFile $DBDIR/$db.nsq temp
	registerFile $DBDIR/$db.nin temp
	registerFile $DBDIR/$db.nhr temp

    fi
	
    blast_tbn=$SEARCHES/${SEEDQUERYFILE}_vs_${db}.tbn.m8
    if needsUpdate $blast_tbn $SEEDQUERYFILE $DBDIR/$db $BLASTBINDIR/blastall
    then
	$BLASTBINDIR/blastall -p tblastn -i $SEEDQUERYFILE -d $DBDIR/$db -m8 -o $blast_tbn;
	registerFile $blast_tbn temp
    fi
    
    blast_tbn_names=$SEARCHES/${SEEDQUERYFILE}_vs_${db}.tbn.e01contigs.names
    if needsUpdate $blast_tbn_names $blast_tbn
    then
	awk '$11 < 0.1 { print; }' < $blast_tbn |cut -f2 |sort |uniq > $blast_tbn_names 
	registerFile $blast_tbn_names
    fi

    blast_tbn_fasta=$SEQTEMP/${SEEDQUERYFILE}_vs_${db}.tbn.e01contigs.fasta
    if needsUpdate $blast_tbn_fasta $blast_tbn_names $DBDIR/${db} $BINDIR/fasta_header_grep.pl
    then
	blast_tbn_names_grepsafe=$SEARCHES/${SEEDQUERYFILE}_vs_${db}.tbn.e01contigs.names.grepsafe
	perl -ne 's/\|$//; s/\|/\\\|/g; s/\./\\\./; print;' < $blast_tbn_names > $blast_tbn_names_grepsafe
	registerFile $blast_tbn_names_grepsafe temp
	$BINDIR/fasta_header_grep.pl -w -f $blast_tbn_names_grepsafe $DBDIR/${db} > $blast_tbn_fasta
	registerFile $blast_tbn_fasta result
    fi

    # split any hits larger than N into n large chunks around hits to simplify life for genewise? Nah, simply truncate all input genomes first..

    wise=$SEARCHES/${MOTIF}_in_${SEEDQUERYFILE}_vs_${db}.tbn.e01contigs.wise
    if needsUpdate $wise $blast_tbn_fasta
    then
	genewisedb $MOTIF $blast_tbn_fasta -aln 300 -cut 10 -pep -sum -hmmer -dnadb -gene worm.gf > $wise
	registerFile $wise result
    fi

    wise_hitsfasta=${wise}.hitsfasta
    if needsUpdate $wise_hitsfasta $wise
    then
	grep -v '^>Results' $wise| grep -v '^Making' |perl -e 'my $fasta = ""; $score_upcoming = 0; $past_alignments_barrier=0; while(my $l = <STDIN>) { chomp $l; if(!$past_alignments_barrier) { if($l=~m/\#Alignments/) {$past_alignments_barrier = 1;} } else { if($fasta ne "" ) { if ($l =~ m/^>(.+)/) { print $fasta; $start=shift @start; $end=shift @end; $score=shift @score; $fasta = $l."_".$start."-".$end."_".$score."\n"; } elsif ($l =~m/\/\//) { print $fasta; $fasta=""; } else { $fasta .= $l."\n"; } } elsif ($l =~m/^>(.+)/) { $start=shift @start; $end=shift @end; $score=shift @score; $fasta = $l."_".$start."-".$end."_".$score."\n" }; if($score_upcoming == 1) { if($l =~m/\/\//) {$score_upcoming=0;} else { @l = split (/\s+/, $l); push @score,$l[0]; push @start,$l[5]; push @end,$l[6]; } } if ($l =~ /^Bits/) { $score_upcoming = 1;} } }' |sed -e 's/.pep//;' > $wise_hitsfasta
	registerFile $wise_hitsfasta temp
    fi

    # clustalnames with a db initial. This is the ugliest thing I've seen in a long time. 

    wise_tree_fasta=$SEARCHES/${MOTIF}_in_${SEEDQUERYFILE}_vs_${db}.tbn.e01contigs.wise.tree.fasta
    if needsUpdate $wise_tree_fasta $wise_hitsfasta
    then
	grep $db $DBDIR/db_initials |cut -f2 > $SEQTEMP/db_init
	registerFile $SEQTEMP/db_init temp

	cat $SEQTEMP/db_init $wise_hitsfasta | perl -e '$l = <STDIN>; chomp $l; $dbinitial=$l; while($l=<STDIN>) { chomp $l; if($l=~/^>[Cc]{1}ontig.*/) {$out=$dbinitial."c"; $l=~s/[Cc]{1}ontig/$out/; } elsif($l=~/^>[Ss]{1}caffold.*/) {$out=$dbinitial."s"; $l=~s/[Ss]{1}caffold/$out/; } elsif($l=~/^>[Ss]{1}upscaffold.*/) {$out=$dbinitial."s"; $l=~s/[Ss]{1}upscaffold/$out/; } elsif($l=~/^>[Ss]{1}uper_.*/) {$out=$dbinitial."s"; $l=~s/[Ss]{1}uper_/$out/; } elsif($l=~/^>[Gg]{1}ene[Ss]{1}caffold.*/) { $out=$dbinitial."gs"; $l=~s/[Gg]{1}ene[Ss]{1}caffold/$out/; } elsif($l=~/^>[Ss]{1}uper[Cc]{1}ontig.*/) {$out=$dbinitial."s"; $l=~s/[Ss]{1}uper[Cc]{1}ontig/$out/; } elsif($l=~/^>Smp_contig.*/) {$out=$dbinitial."c"; $l=~s/Smp_contig/$out/; } elsif($l=~/^>Bmal_supercontig.*/) {$out=$dbinitial."s"; $l=~s/Bmal_supercontig/$out/; } elsif($l=~/^>supercont/) { $out=$dbinitial."s"; $l=~s/>supercont/>$out/; } elsif($l=~/^>chr/) { $out=$dbinitial."ch"; $l=~s/>chr/>$out/; } elsif($l=~/^>[Gg]{1}roup/) { $out=$dbinitial."g"; $l=~s/>[gG]roup/>$out/; } elsif($l=~/(\d+)_countedlength_\d+/ ) { $out=$dbinitial."_".$1; $l=~s/\d+_countedlength_\d+/$out/} elsif($l=~m/gi\|\d+\|(?:emb|gb)\|/) { $out=$dbinitial."_"; $l=~s/gi\|\d+\|(?:emb|gb)\|/$out/; $l=~s/\|_/_/g;} elsif($l=~/^>/) { $l=~s/^>/>$dbinitial/; $l=~s/CAAB010//g; $l=~s/GJ0//g; $l=~s/_unplaced//g; }; print $l,"\n";}' > $wise_tree_fasta
    fi
done

if [ "$DOTREES" -eq 1 ]
then
#complete tree
    cp ${SEQTEMP}/${MOTIF}_in_${SEEDQUERYFILE}.clustalnames $SEARCHES/${MOTIF}_in_${SEEDQUERYFILE}.tbn.e01contigs.wise.complete.tree.fasta
    for db in $DBS ; do 
	cat $SEARCHES/${MOTIF}_in_${SEEDQUERYFILE}_vs_${db}.tbn.e01contigs.wise.tree.fasta >> $SEARCHES/${MOTIF}_in_${SEEDQUERYFILE}.tbn.e01contigs.wise.complete.tree.fasta
    done
    clustalw2 -INFILE=${SEARCHES}/${MOTIF}_in_${SEEDQUERYFILE}.tbn.e01contigs.wise.complete.tree.fasta
    clustalw2 -INFILE=${SEARCHES}/${MOTIF}_in_${SEEDQUERYFILE}.tbn.e01contigs.wise.complete.tree.aln -BOOTSTRAP=1000

# per db trees
#    for db in $DBS ; do 
#	cat ${SEQTEMP}/${MOTIF}_in_${SEEDQUERYFILE}.clustalnames >> $SEARCHES/${MOTIF}_in_${SEEDQUERYFILE}_vs_${db}.tbn.e01contigs.wise.tree.fasta
#	clustalw $SEARCHES/${MOTIF}_in_${SEEDQUERYFILE}_vs_${db}.tbn.e01contigs.wise.tree.fasta
#	clustalw -INFILE=$SEARCHES/${MOTIF}_in_${SEEDQUERYFILE}_vs_${db}.tbn.e01contigs.wise.tree.aln -BOOTSTRAP=1000
#   done
fi

#done for now..

exit;

# subtree select sequences once
rm $SEARCHES/${MOTIF}_in_${SEEDQUERYFILE}.tbn.e01contigs.wise.phb.select.acr23.dbnames.fasta

for db in $DBS ; do 
    ./get_subtrees.pl $SEARCHES/${MOTIF}_in_${SEEDQUERYFILE}_vs_${db}.tbn.e01contigs.wise.tree.phb selectlist_acr23 $SEARCHES/${MOTIF}_in_${SEEDQUERYFILE}_vs_${db}.tbn.e01contigs.wise.tree.phb.selectlist_acr23
    perl -ne 's/\|/\\\|/g; print;' < $SEARCHES/${MOTIF}_in_${SEEDQUERYFILE}_vs_${db}.tbn.e01contigs.wise.tree.phb.selectlist_acr23 >  $SEARCHES/${MOTIF}_in_${SEEDQUERYFILE}_vs_${db}.tbn.e01contigs.wise.tree.phb.selectlist_acr23.grepsafe
    ./fasta_header_grep.pl -f $SEARCHES/${MOTIF}_in_${SEEDQUERYFILE}_vs_${db}.tbn.e01contigs.wise.tree.phb.selectlist_acr23.grepsafe $SEARCHES/${MOTIF}_in_${SEEDQUERYFILE}_vs_${db}.tbn.e01contigs.wise.tree.fasta > $SEARCHES/${MOTIF}_in_${SEEDQUERYFILE}_vs_${db}.tbn.e01contigs.wise.tree.phb.selectlist_acr23.fasta
    #unique initial  --moved upstream
done

perl -ne 'chomp; if (/^>(\S+)/) {print $1,"\n";}' $SEARCHES/${MOTIF}_in_${SEEDQUERYFILE}.tbn.e01contigs.wise.phb.select.acr23.dbnames.fasta |sort|uniq >$SEARCHES/${MOTIF}_in_${SEEDQUERYFILE}.tbn.e01contigs.wise.phb.select.acr23.uniqnames

perl -ne 's/\|/\\\|/g; print;' < $SEARCHES/${MOTIF}_in_${SEEDQUERYFILE}.tbn.e01contigs.wise.phb.select.acr23.uniqnames > $SEARCHES/${MOTIF}_in_${SEEDQUERYFILE}.tbn.e01contigs.wise.phb.select.acr23.uniqnames.grepsafe

./fasta_header_grep.pl -f $SEARCHES/${MOTIF}_in_${SEEDQUERYFILE}.tbn.e01contigs.wise.phb.select.acr23.uniqnames.grepsafe $SEARCHES/${MOTIF}_in_${SEEDQUERYFILE}.tbn.e01contigs.wise.phb.select.acr23.dbnames.fasta > $SEARCHES/${MOTIF}_in_${SEEDQUERYFILE}.tbn.e01contigs.wise.phb.select.acr23.uniq.fasta

# perl -ne 's/gi\|\d+\|//g; s/\|_/_/g; print; ' <  $SEARCHES/${MOTIF}_in_${SEEDQUERYFILE}.tbn.e01contigs.wise.phb.select.acr23.uniq.fasta > $SEARCHES/${MOTIF}_in_${SEEDQUERYFILE}.tbn.e01contigs.wise.phb.select.acr23.uniq.clustalnames

clustalw $SEARCHES/${MOTIF}_in_${SEEDQUERYFILE}.tbn.e01contigs.wise.phb.select.acr23.uniq.fasta
clustalw -INFILE=$SEARCHES/${MOTIF}_in_${SEEDQUERYFILE}.tbn.e01contigs.wise.phb.select.acr23.uniq.aln -BOOTSTRAP=1000

#oops on CJA?! others? acr-5?! nooooot complete yet..

#for db in $DBS ; do for name in `cat ${SEEDQUERYFILE}_vs_${db}.tbn.e001contigs.names`; do seqret ../${db}:${name} ../seq_temp/${db}.${name}.fasta ; done ; done
#for db in $DBS ; do rm ${SEEDQUERYFILE}_vs_${db}.tbn.e001contigs.fasta; for contig in `cat searches/${SEEDQUERYFILE}_vs_${db}.tbn.e001contigs.names`; do cat seq_temp/${db}.${contig}.fasta >> seq_temp/${SEEDQUERYFILE}_vs_${db}.tbn.e001contigs.fasta; done ; done
#extremely slow. try an @file USA or use own fastaheadergrep instead?

# FOR -trans
#    grep -v '^>Results' searches/${MOTIF}_in_${SEEDQUERYFILE}_vs_${db}.tbn.e001contigs.wise | perl -e 'my $fasta = ""; $score_upcoming = 0; while(my $l = <STDIN>) { chomp $l; if($fasta ne "" ) { if ($l =~ m/^>(.+)/) { print $fasta; $fasta = $l."\n"; } elsif ($l =~m/\/\//) { print $fasta; $fasta=""; } else { $fasta .= $l."\n"; } } elsif ($l =~m/^>(.+)/) { $fasta = $l."\n" }; }'|sed -e 's/.sp.tr//; s/\[//; s/\]//; s/:/-/' > searches/${MOTIF}_in_${SEEDQUERYFILE}_vs_${db}.tbn.e001contigs.wise.hitsfasta ;

#    genewisedb Neur_chan_LBD_ls.hmm seq_temp/${SEEDQUERYFILE}_vs_${db}.tbn.e001contigs.fasta -trans -sum -hmmer -dnadb -gene worm.gf > searches/${MOTIF}_in_${SEEDQUERYFILE}_vs_${db}.tbn.e001contigs.wise ; 

#    grep -v -f selectlist_acr23 searches/${MOTIF}_in_${SEEDQUERYFILE}_vs_${db}.tbn.e01contigs.wise.tree.phb.selectlist_acr23 > searches/${MOTIF}_in_${SEEDQUERYFILE}_vs_${db}.tbn.e01contigs.wise.tree.phb.selectlist_acr23_noce
#./fasta_header_grep.pl -f selectlist_acr23 Neur_chan_LBD_ls.hmm_vs_ce_go0005230.clustalnames > searches/${MOTIF}_in_${SEEDQUERYFILE}.tbn.e01contigs.wise.phb.select.acr23.fasta

#       cat searches/${MOTIF}_in_${SEEDQUERYFILE}_vs_${db}.tbn.e01contigs.wise.hitsfasta Neur_chan_LBD_ls.hmm_vs_ce_go0005230.clustalnames |perl -ne 's/gi\|\d+\|//g; s/\|_/_/g; print; ' > searches/${MOTIF}_in_${SEEDQUERYFILE}_vs_${db}.tbn.e01contigs.wise.tree.fasta

#    clustalw searches/${MOTIF}_in_${SEEDQUERYFILE}_vs_${db}.tbn.e01contigs.wise.tree.clustalnames

#    grep $db db_initials |cut -f2 > db_init
#    cat db_init searches/${MOTIF}_in_${SEEDQUERYFILE}_vs_${db}.tbn.e01contigs.wise.tree.phb.selectlist_acr23.fasta | perl -e '$l = <STDIN>; chomp $l; $dbinitial=$l; while($l=<STDIN>) { chomp $l; if($l=~/^>.*([Cc]{1}ontig\S+).*/) {$out=$dbinitial."_Contig"; $l=~s/[Cc]{1}ontig/$out/; } elsif($l=~/^>chr/) { $out=$dbinitial."_chr"; $l=~s/>chr/>$out/; } elsif(m/gi\|\d+\|emb\|/) { $out=$dbinitial; $l=~s/gi\|\d\+\|emb\|/>$out/; }; print $l,"\n";}' >> searches/${MOTIF}_in_${SEEDQUERYFILE}.tbn.e01contigs.wise.phb.select.acr23.dbnames.fasta
