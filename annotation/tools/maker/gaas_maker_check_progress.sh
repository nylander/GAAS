#!/bin/bash

############################################################################
# JD 2014/04
# What it does ? => script look the index log of maker to resume the contig analyzed and those not
# use:  script logfile
############################################################################


# Arguments and Paths
logFile=""
if (( $# !=1 )); then
	echo -e "You can specify the maker log file path.\nNevetheless we will try the default path ..."
	if [ ! -f genome.maker.output/genome_master_datastore_index.log ];then
		echo -e "The default path do not work\nRetry giving a correct path to the log file"
		exit		
	else
		logFile=genome.maker.output/genome_master_datastore_index.log
	fi
else
	logFile=$1	
fi


# Light information
#small=$(awk ' { if($3 == "SKIPPED_SMALL") nb++} END {print nb}' $logFile)
small=$(awk ' { if($3 == "SKIPPED_SMALL") print $1}' $logFile | sort -n | uniq | wc -l)
if [[ $small == "" ]];then
	small=0
fi
echo -e "\n$small contigs are too small to be analyzed\n"

#startedDoublon=$(awk ' { if($3 == "STARTED") nb++} END {print nb}' $logFile)
startedOne=$(awk ' { if($3 == "STARTED") print $1}' $logFile | sort -n | uniq -u | wc -l)
echo "$startedOne contigs has begin to be studied only one time"
startedDoublon=$(awk ' { if($3 == "STARTED") print $1}' $logFile | sort -n | uniq -d | wc -l)
echo "$startedDoublon contigs have been started several times"
((started=$startedOne+$startedDoublon))
echo "$started contigs has begin to be studied in total"
startedAll=$(awk ' { if($3 == "STARTED") print $1}' $logFile | wc -l)
echo -e "$startedAll STARTED signal in total\n"

#finished=$(awk ' { if($3 == "FINISHED") nb++} END {print nb}' $logFile)
finished=$(awk ' { if($3 == "FINISHED") print $1}' $logFile | sort -n | uniq -u | wc -l)
finishedDoublon=$(awk ' { if($3 == "FINISHED") print $1}' $logFile | sort -n | uniq -d | wc -l)
echo "$finished contigs have been finished"
echo -e "$finishedDoublon contigs have been finished several times\n"

#skippedperm=$(awk ' { if($3 == "DIED_SKIPPED_PERMANENT") nb++} END {print nb}' $logFile)
skippedpermAll=$(awk ' { if( $3 ~ /.*DIED.*/ ) print $1}' $logFile | wc -l)
echo "$skippedpermAll DIED_SKIPPED_PERMANENT signal found. The number of retry attempts has been reached"
skippedperm=$(awk ' { if( $3 ~ /.*DIED.*/ ) print $1}' $logFile | sort -n | uniq | wc -l)
echo -e "$skippedperm contigs (doublon removed) have the signal DIED_SKIPPED_PERMANENT. The number of retry attempts has been reached\n"

declare -A startList
declare -A finishList
declare -A diedList
declare -A allBug
for i in $(awk ' { if($3 == "STARTED") print $1}' $logFile | sort -n | uniq);do
	startList["$i"]=1
done

for i in $(awk ' { if($3 == "FINISHED") print $1}' $logFile | sort -n | uniq);do
	finishList["$i"]=1
done

for i in $(awk ' { if($3 ~ /.*DIED.*/) print $1}' $logFile | sort -n | uniq);do
        diedList["$i"]=1
	allBug["$i"]=1
done

# Verif	STARTED signal whithout FINISHED signal
cptNeverFinish=0
for i in ${!startList[*]};do
        if [[ ! ${finishList[$i]-X} == ${finishList[$i]} ]];then
#                echo "$i never finished ..."
#                grep "$i/" $logFile
                ((cptNeverFinish=cptNeverFinish+1))
        	allBug["$i"]=1
	fi
done
echo -e "$cptNeverFinish contig whithout FINISHED signal.\n"


# Verif FINISHED signal whithout STARTED signal
cptNeverStart=0
for i in ${!finishList[*]};do
        if [[ ! ${startList[$i]-X} == ${startList[$i]} ]];then
#		echo "$i never started ..."
# 		grep "$i/" $logFile
		((cptNeverStart=cptNeverStart+1))
        fi
done
echo -e "$cptNeverStart contig whithout STARTED signal.\n"

#About All errors detected
echo -e "\n###Finally, ${#allBug[@]} errors has been found in the log file:(See below)###"
num=1
for i in ${!allBug[*]};do
	echo -e "\n$num) $i:"
	grep "$i/" $logFile
	((num=num+1))
done

# verify in directory the errors
bugcounter=${#allBug[@]}
echo -e "\n### Now, verification of errors found in the directory of these ananlysis ###"
for i in ${!allBug[*]};do
        echo -ne "\n$i:"
	IFS=$'\n'
        for j in $(grep "$i/" $logFile);do #for each line containing the contig bugged
		if [[ $j =~ [^[:alnum:]+[:blank:].+[:blank:].+]  ]];then #if the line is correctly written. Mean should have, at least, three columns
			pathdir=$(echo $j | awk '{print $2}') # get path in the second column
			goodPath="$(dirname $logFile)/$pathdir/run.log" # get the root path and concatenate to the path took before 
			ok="no"
			for k in $(grep "FINISHED" $goodPath);do #get all line with FINISHED term
				if [[ $(echo $k | grep "final.section") != "" ]];then # verify if one of these line contain final.section
					echo -ne "Finaly this analysis was well terminated !"
					ok="yes"
					((bugcounter=bugcounter-1)) #We remove this case of the counter because it is not a bug
					break
				fi
			done
			if [[ $ok == "no" ]];then #if no line contains Finished or final.section 
# Print Pretty
				sizeString=$(echo $i | wc -m)
				lengthBeforeWritte=25;
				nbSpaceToWrite=$((lengthBeforeWritte-sizeString))
				startValue=1
				while [[ $startValue -le $nbSpaceToWrite ]];do
					echo -ne "-";((startValue=startValue+1))
				done
				echo -ne "It seems that this analysis has never finished"
				if [[ $(grep "DIED" $goodPath) != "" ]];then
					echo "Indeed the DIED signal is present (see below)"
					grep "DIED" $goodPath
				fi

			fi
			break
		fi
	done
done

if (( $bugcounter == 0 ));then
	echo -e "\n\nThis job does not contains bug ! Congratulation."
else
	echo -e "\n\nWe found $bugcounter errors. Good luck for next ..."
fi

echo -e "\nBye"
