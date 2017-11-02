#!/bin/bash

#This script
cleanIntermediateFile="yes"

if [[ $1 == "no" ]];then
	cleanIntermediateFile="no"
fi

for i in {0..50}_*;do

	if [[ -f $i ]];then

		if [[ ! $i =~ ^[[:digit:]]+_correct ]];then
			echo -e "\nTest of $i";
			testperfect="no"
			~/git/NBIS/GAAS/annotation/Tools/Converter/gxf_to_gff3.pl --gff $i -o test.gff3  &> /dev/null  
			pref=$(echo $i | cut -d'_' -f1)
			fileok=${pref}_correct_output.gff

			if [ ! -f $fileok ];then
				echo "We didnt find any correct output to check against for $i ( $fileok ) "
			else
				resu=$(diff test.gff3 $fileok)
				if [[ $resu != "" ]];then
					echo -e "There is differences between the correct reference output and the current output:\n$resu"
				else
					echo "test1 ok !"
					testperfect="yes"
				fi
			fi

			#echo "check against itself"
			~/git/NBIS/GAAS/annotation/Tools/Converter/gxf_to_gff3.pl --gff test.gff3 -o test2.gff3  &> /dev/null  
			resu=$(diff test2.gff3 $fileok)
			if [[ $resu != "" ]];then
					echo -e "There is differences between the original current output and the output of this file processed again:\n$resu"
			else
					echo "test2 ok !"
					if [[ $testperfect == "yes" ]];then
						echo "All test perfect !"
					fi
			fi
			
		fi
	fi
done

#Intermediate file cleaned
if [[ $cleanIntermediateFile == "yes" ]];then
	rm test.gff3
	rm test2.gff3
fi

