#!/bin/bash

############################################################################
# JD 2014/04
# What it does ? read gff file to define te number of gene and mRNA and the AED mean
# use:  script gffFile 
############################################################################

# Arguments and Paths
if (( $# !=1 )); then
        echo -e "You have to specify a a gff file"
	exit
else
	gffFile=$1
	if [[ ! -f $gffFile ]];then 
		echo "The file $gffFile doesn't exist."
		exit
	fi
fi

#file arg is required
total=0;nb=0; for i in $(awk '{if($3 == "mRNA") print $9}' $1 | cut -d';' -f4 | cut -d'=' -f2 );do resu=$total; total=$(echo $resu+$i | bc);((nb=nb+1));done; mean=$(echo "scale=4;$total / $nb" | bc)
nbGene=$(awk '{if($3 == "gene") print $0}' $1 | wc -l)
echo -e "\nThere is $nbGene genes for $nb mRNA."
echo -e "\nAED=> Between 1 and 0. The lowest the this value is, the better it is."
echo "AED Total of $nb mRNA is $total"
echo "The AED mean by mRNA is $mean"
