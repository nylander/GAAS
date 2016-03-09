#!/bin/bash

############################################################################
# JD 2014/04
# What it does ? => script look the index log of maker to resume the contig analyzed and those not
# use:  script logfile
############################################################################


# Arguments and Paths
dir0=""
if (( $# !=1 )); then
	echo -e "You can specify the maker datastore path.\nNevetheless we will try the default path ..."
	if [ ! -f genome.maker.output/genome_datastore ];then
		echo -e "The default path do not work\nRetry giving a correct path to the datastore directory"
		exit		
	else
		dir0=genome.maker.output/genome_datastore
	fi
else
	dir0=$1	
fi

##Function
# verify in errors in directory 

function testContigState {
        pathContig=$1
	goodPath=${pathContig}run.log
	ok="no"
	if [[ ! $(tail -n1 $goodPath | grep "trnascan" ) == "" ]];then
#		echo "nothing annotated in this one"
		((skippedSmall=skippedSmall+1))
	else {
#		if [[ $(grep "log.child" $goodPath) != "" ]];then
#			echo "This studies have child.log $pathContig"
#		fi	

	        for k in $(grep "FINISHED" $goodPath);do #get all line with FINISHED term
        		if [[ $(echo $k | grep "final.section") != "" ]];then # verify if one of these line contain final.section
                        	ok="yes"
                        	((analyzeOk=analyzeOk+1))
                        	break
                	fi
        	done
        	if [[ $ok == "no" ]];then #if no line contains Finished or final.section
                        echo "It seems that this analysis has never finished $pathContig"
			if [[ $(grep "DIED" $goodPath) != "" ]];then
                        	echo "Indeed the DIED signal is present (see below)"
                                grep "DIED" $goodPath
				((skippedDied=skippedDied+1))
                        else {
				echo "/!\\ Wierd case. You have to verify manually if this study was perform correctly."
			}
			fi
                fi
	}
	fi
}

# Light information
nbRlevel1=$(ls -l -d ${dir0}*/ |wc -l)
echo "There is $nbRlevel1 directory at level 1"

skippedSmall=0
analyzeOk=0
skippedDied=0

counterL2=0
counterL3=0
#counterL3cl=0
counterT=0
IFS=$'\n'
max=100000

for dirLevel2 in $(ls -l -d ${dir0}*/ | awk '{ print $9 }' );do
	if [[ $dirLevel2 == "" ]];then
		continue
	fi

	nbRlevel2=$(ls -l -d ${dirLevel2}*/ | wc -l)
	((counterL2=counterL2+nbRlevel2))
	((counterT=counterT+nbRlevel2))
	echo "There is $nbRlevel2 directory at level 2 for $dirLevel2"

	for dirLevel3 in $(ls -l -d ${dirLevel2}*/ | awk '{ print $9 }' );do
		 if [[ $dirLevel3 == "" ]];then
	                continue
       		 fi
		
        	nbRlevel3=$(ls -l -d $dirLevel3*/ | wc -l)
		((counterL3=counterL3+nbRlevel3))		
		for dirFinalLevel in $(ls -l -d ${dirLevel3}*/ | awk '{ print $9 }' );do
			testContigState $dirFinalLevel	
		done
#		if (($nbRlevel3 != 1 ));then
#			echo "    =>There is $nbRlevel3 directory at level3 for $dirLevel3"
#			((counterL3cl=counterL3cl+nbRlevel3-1))
#			((counterT=counterT+nbRlevel3-1))
#		fi
		if (( $counterL3 >= $max ));then
			break
		fi
	
	done
	
	if (( $counterL3 > $max ));then
                break
        fi
	
done

echo -e "\n#Resume:"
echo -e "In total, we have $nbRlevel1 directory at level 1"
echo -e "In total, we have $counterL2 directory at level 2"
#echo -e "In total, we have $counterL3cl directory at level 3 clean"
echo -e "In total, we have $counterL3 directory at level 3 (correspond to thenimber of Contig studied)"
#echo -e "In total, we have $counterT directory in total"


#Resume result // Display them
echo -e "\nAmong $counterL3 Contigs studied we have:"
echo "$skippedSmall contig(s) too small"
echo "$skippedDied contig(s) not studied (Number attemps reaching)"
echo "$analyzeOk contig(s) studied succesfully"

((totalVerified=skippedSmall+analyzeOk+skippedDied))
if (( $totalVerified != $counterL3 ));then
	((missingInfo=$counterL3-$totalVerified))
	echo "There is a problem. We didn't understand what happened for $missingInfo"
fi 

echo -e "\nBye"
