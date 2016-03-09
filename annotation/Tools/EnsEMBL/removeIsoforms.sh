#!/bin/bash

# JD 2015
#
#This script kept only one isoform per gene for proteomes coming form Ensembl
####

# Arguments and Paths
############################################################################

if (( $# != 2 )) ; then
	echo -e "The script allows to filter proteome in fasta format from Ensembl with aims to keep the longest isoform per gene !"
        echo -e "The script needs 2 parameters: \n(1)The proteome input fasta file"
        echo -e "(2)The cleaned proteome output in fasta format"
        exit
fi


 cat $1 | awk '/^>/ {if(N>0) printf("\n"); printf("%s\t",$1"\t"$4);N++;next;}{printf("%s",$0);} END {if(N>0) printf("\n");}' | awk -F '\t'  '{printf("%s\t%d\n",$0,length($3));}' | sort -t '	' -k2,2 -k4,4nr  | sort -t '	' -k2,2 -u -s | cut -f 1,2,3 | awk '{print $1"\t"$2"\n"$3}' | fold -w 60 > $2
	








# Explanation
#############
# cat test | awk '/^>/ {if(N>0) printf("\n"); printf("%s\t",$1"\t"$4);N++;next;}{printf("%s",$0);} END {if(N>0) printf("\n");}' |\ #linearize fasta and print col1 col4 and seq linearized
# awk -F '\t'  '{printf("%s\t%d\n",$0,length($3));}' |\ #extact length on the 4th column
# sort -t '      ' -k2,2 -k4,4nr  |\#sort on column2, inverse length
# sort -t '    ' -k2,2 -u -s |\# #sort on column 2, unique, stable sort (keep previous order)
# cut -f 1,2,3 |\ #cut 3column
# awk '{print $1"\t"$2"\n"$3}' |\ # print col1\tcol2 \n col3
# fold -w 60 #pretty fasta
