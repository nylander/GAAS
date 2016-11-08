#!/bin/bash

for i in {1..20}_*;do
	if [[ -f $i ]];then
		if [[ ! $i =~ ^[1-9]+_correct ]];then
			echo -e "\nTest of $i";
			~/git/NBIS/GAAS/annotation/Tools/Util/gff/gff3_sp_IO.pl --gff $i -o test.gff  &> /dev/null  
			pref=$(echo $i | cut -d'_' -f1)
			fileok=${pref}_correct_output.gff

			if [ ! -f $fileok ];then
				echo "We didnt find any output to check against for $i ( $fileok ) "
			else
				resu=$(diff test.gff $fileok)
				if [[ $resu != "" ]];then
					echo "$resu"
				else
					echo "Result perfect !"
				fi
			fi

			echo "check against itsefl"
			~/git/NBIS/GAAS/annotation/Tools/Util/gff/gff3_sp_IO.pl --gff test.gff -o test2.gff  &> /dev/null  
			resu=$(diff test2.gff $fileok)
			if [[ $resu != "" ]];then
					echo "$resu"
			else
					echo "Result perfect !"
			fi
		fi
	fi
done
rm test.gff
rm test2.gff
