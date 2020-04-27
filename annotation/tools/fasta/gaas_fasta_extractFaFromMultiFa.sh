#!/bin/bash

# 10/2014
# Jacques Dainat
# jacques.dainat@nbis.se
#

#ARGS in command line
if (( $# < 2 )); then
        echo -e "This script allows to extract easily and efficiently the sequence(s) from multi-fasta-file\n"
	echo -e "The script needs 2 parameters"
        echo -e "[usage: script.sh multiFastaFile.fa [HeaderName or HeaderNameFile] [m f] ]"
        echo -e "In the case of HeaderNameFile, this file must contains 1 header per line"
        echo -e "m <= term match\nf <= line finish by this term"
	exit
else
        multiFasta=$1
        name=$2
	matchType=$3
fi

if [[ $matchType == "" ]];then
	echo -e "By default we are looking for matching term"
	matchType="m"
fi

#PREPARE REGEX
regex="$name";
if [[ $matchType == "m" ]];then
	regex="$name"
elif [[ $matchType == "f" ]];then
	regex="${name}\$"
fi

if [[ ! -e $name ]];then
	nbOcc=$(grep -c "$regex" $multiFasta)
	if [[ $nbOcc == 0 ]];then
                echo "No match found"; exit 1;
        elif [[ $nbOcc > 1 ]];then
		echo "We found $nbOcc match. Only one is mandatory"; exit 1;
	fi

    line=$(grep -nr "$regex" $multiFasta)
    lineNumber=$(echo $line | cut -d':' -f1)
    sed -n "$lineNumber p" $multiFasta;
    awk -v nb=$lineNumber 'NR > nb {if ($0 ~ ">") exit; else print $0 }' $multiFasta
## Case of header file
else
    echo "Reading $name file that should contain one header name by line."
    IFS=$'\n';
    for i in $(cat $name);do
        line=$(grep -nr "${i}" $multiFasta);
            if [[ ! -z $line ]];then
                for j in $line;do
                    lineNumber=$(echo $j | cut -d':' -f1)
                    sed -n "$lineNumber p" $multiFasta
                    awk -v nb=$lineNumber 'NR > nb {if ($0 ~ ">") exit; else print $0 }' $multiFasta
                done
            fi
    done
fi

exit

