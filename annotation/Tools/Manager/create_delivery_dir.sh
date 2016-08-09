#!/bin/bash

# JD 2015
#
#This script will soft-mask (lowerCase) a genome using a gff file as input
####



# Arguments and Paths
############################################################################
dirname=""

if (( $# != 1 )) ; then
        echo -e "No directory name given. By default we will call it <delivery> "
	dirname="delivery"
else
	if [[ $1 =~ \-?[hH]{1}elp$ || $1 =~ \-?[hH]{1}$ ]];then
		echo -e "This script prepre the tree of directories for genome annotation results as well as a Readme file.\nCommand:\n========\nscript.sh [dirName]\nBy default the directory name is called <delivery>"
		exit
	else	
		dirname=$1
	fi
fi

if [[ -d $dirname ]];then
	echo "file <$dirname> already exists"; exit;
else
	mkdir -p "$dirname/metadata"
	mkdir -p "$dirname/summary"
	mkdir -p "$dirname/gff/repeats"
	mkdir -p "$dirname/gff/lift-over"
	mkdir -p "$dirname/gff/ncRNA"
	mkdir -p "$dirname/gff/gene-build"
	mkdir -p "$dirname/gff/transcript"
	mkdir -p "$dirname/fasta/transcript"
	
	#write Readme
	echo -e "* BILS Genome Annotation Platform team *" >> $dirname/ReadMe.txt
	mydate=$(date)
	echo -e "$mydate\nPlease find here an overview about data available.\n" >> $dirname/ReadMe.txt
	echo -e "<Species_Annotation_Report.docx> file => A complete report of work/process performed." >> $dirname/ReadMe.txt
	echo -e "================================\n" >> $dirname/ReadMe.txt
	echo -e "<fasta> directory => It contains data in FASTA format." >> $dirname/ReadMe.txt
	echo -e "=================	<cds.fasta> file => It contains cds of annotated genes." >> $dirname/ReadMe.txt
	echo -e "			<protein.fasta> file => It contains the proteins produced by the annotated genes." >> $dirname/ReadMe.txt
	echo -e "	<transcripts> sub directory" >> $dirname/ReadMe.txt
	echo -e "	===========     <Trinity_*_.fasta> is the output of our Trinity pipeline used as inout for Maker.\n" >> $dirname/ReadMe.txt
	echo -e "<gff> directory => It contains data in gff3 format." >> $dirname/ReadMe.txt
	echo -e "===============\n" >> $dirname/ReadMe.txt
	echo -e "	<gene-builds> sub directory\n" >> $dirname/ReadMe.txt
	echo -e "	===========	<species_rcX.gff> Are the “release candidate” annotation for the species. The higher release candidate corresponds to the most successful achievement.\n" >> $dirname/ReadMe.txt
	echo -e "	<lift-overs> sub directory (depending on customer demand)" >> $dirname/ReadMe.txt
	echo -e "	============	<referenceSpecies2studiedSpecies.gff> It is the lift-over of a reference species genome annotation (from Ensembl) on the genome under investigation.\n" >> $dirname/ReadMe.txt
	echo -e "	<ncRNA> sub directory" >> $dirname/ReadMe.txt
	echo -e "	=======		<ncRNA_rfam.gff> It contains annotation of non-coding element annotated using the Eucaryotes data from Rfam database." >> $dirname/ReadMe.txt
	echo -e "			<*_tRNA.gff> It contains transfert RNA annotated using tRNAscan.\n" >> $dirname/ReadMe.txt
	echo -e "	<repeats> sub directory" >> $dirname/ReadMe.txt
	echo -e "	=======	<repeatmasker.gff> It contains repeats annotated thanks to repeatmasker." >> $dirname/ReadMe.txt
	echo -e "		<repeatrunner.gff> It contains repeats annotated thanks to repeatrunner.\n" >> $dirname/ReadMe.txt
	echo -e "	<transcripts> sub directory" >> $dirname/ReadMe.txt
	echo -e "	===========       <Cufflinks_*.gff> contains output of our Cufflinks pipeline used as input for Maker. (Depending of RNAseq data types, and the most appropriate Transcriptome assembly method chosen.)" >> $dirname/ReadMe.txt
	echo -e "			  <Trinity_*.gff> is the gff output created by Maker using the Trinity_*_.fasta file as input. (Depending of RNAseq data types, and the most appropriate Transcriptome assembly method chosen.)" >> $dirname/ReadMe.txt
	echo -e "			  <*> Depending of RNAseq data types, and the most appropriate Transcriptome assembly method chosen.\n" >> $dirname/ReadMe.txt
	echo -e "<metadata> directory" >> $dirname/ReadMe.txt
	echo -e "==========	<go.txt> Tab-delimited format file (2 rows) containing GO terms of functions retrieved for each transcript." >> $dirname/ReadMe.txt
	echo -e "		<pfam.txt> Tab-delimited format file (2 rows) containing pram terms of functions retrieved for each transcript." >> $dirname/ReadMe.txt
	echo -e "		<interpro.txt> Tab-delimited format file (2 rows) containing interpro terms of functions retrieved for each transcript.\n" >> $dirname/ReadMe.txt
	echo -e "<summary> directory" >> $dirname/ReadMe.txt
	echo -e "=========	<coding_gene.txt> General information about annotated coding genes" >> $dirname/ReadMe.txt
	echo -e "		<tRNA.txt> General information about annotated tRNA" >> $dirname/ReadMe.txt

	echo -e "Well done ! Directories and ReadMe have been created. Think to check The ReadMe..."
	
fi

