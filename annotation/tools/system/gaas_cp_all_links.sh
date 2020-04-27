#!/bin/bash

# JD 2014
#
#This script copy all the links of a directory and create the same link in the selected directory
####



# Arguments and Paths
############################################################################

if (( $# != 2 )) ; then
	echo -e "The script needs 2 parameters: \n(1)Directory where all the link will be copied"
	echo -e "(2)Directory where all the link will be pasted"
	exit
fi 

DirPath=${1%/*}
if [[ $DirPath =~ [\/]$ ]];then
	DirPath=$(echo ${DirPath::-1})
fi

IFS=$'\n'
for i in $(ls -l $1);do 
	if [[ $i =~ ^[l] ]];then
		paired=$(echo $i | awk '{print $9" "$11}') 
		to=$(echo $paired | cut -d' ' -f1)
		from=$(echo $paired | cut -d' ' -f2)
                if [[ $from =~ ^[\/] ]];then
			 ln -s ${from} $2/$to
		else
			ln -s ${DirPath}/${from} $2/$to
		fi
	fi
done 
