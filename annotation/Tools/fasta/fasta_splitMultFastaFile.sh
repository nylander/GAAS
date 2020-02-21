#!/bin/bash

############################################################################
# JD 2014/03
# Cut a multi FastaFile in one file by fasta
# use: script fileInput DirectoryOutpßut 
#############################################################################

#constants
java=$(which java)

# Arguments and Paths
if (( $# !=2 )); then
	echo -e "The script needs 2 parameters: \n(1)The Multi-fasta file as input"
	echo -e "(2)The second corresponds to directory of fasta files output"
	exit
fi 

pathDir=$(pwd)
echo "The directory path is $pathDir"

resuDir=$2
pathresuDir="$pathDir/$resuDir"
[[ $pathresuDir != */ ]] && pathresuDir="$pathresuDir"/ #test if a slash exist at the end
if [ ! -d "$pathresuDir" ]; then
  mkdir $pathresuDir
else
	echo "The directory $pathresuDir already exists !"
	rm -r $pathresuDir
	mkdir $pathresuDir
fi

fastaFile=$1
pathFastaFile="$pathDir/$fastaFile"

awk -v path=$pathresuDir '/^>/{head=gsub(">","");f=$head".fa"; header="true"} {if (header == "true") {print ">"$0 > path"/"f ; header="false" ;} else {print > path"/"f}}' $pathFastaFile
