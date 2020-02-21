#!/bin/bash

## jacques.dainat@nbis.se - 2015
# This script allows to make weekly/monthly/yearly compressed backup.
# It is not completely generalized but can be done easily with some modification.
# In the folder where you want to launch the backup you need a current_backup.log file (last modified by the last backup) and a folder called current_backup containing everything you want to backup

now=$(date +"%Y-%m-%d")
gap=43200 #equivalent to 12 hours. We need it because rsync and backup script are not launched at the same time in order to check the current rsync made bedore is the one done in time. 

# Move to the working directory. Folder where all the backup process and save has to be done
cd /databases/backup/backup_prod/

###############
#### Manage directories needed
weeklyDir="weekly"
monthlyDir="monthly"
yearlyDir="yearly"

if [ ! -d "$weeklyDir" ]; then
	mkdir $weeklyDir
fi
if [ ! -d "$monthlyDir" ]; then
        mkdir $monthlyDir
fi
if [ ! -d "$yearlyDir" ]; then
        mkdir $yearlyDir
fi

#####
# INFO
dayInterval=7
nbWeeklyMax=4
monthInterval=1
nbMonthlyMax=11
yearInterval=1

dateFileCurrentBackup=`date -r current_backup.log  +"%s"`


########
# METHOD

prepare_backup_to_store(){
	
	# copy backup
        mkdir ${now}.backup
        cd ${now}.backup
        ln -s ../current_backup.log
        ln -s ../current_backup
        cd ..
        #compress it
        tar -chzf ${now}.tar.gz ${now}.backup && rm -Rf ${now}.backup
	if [[ $? != 0 ]];then date > READMEplease.log ; echo -e "The exit status code of the tar command is: $? \nPlease check why is not 0! It could be unsafe for the backup.\n" > READMEplease.log;fi
}

####################
#lastBackup // Handle compare to weekly backup
if [ ! "$(ls -A $weeklyDir)" ]; then # directory empty
	prepare_backup_to_store
	# Move it in weekly directory
	mv ${now}.tar.gz $weeklyDir/
else # directory not empty
	nameFileLastBackup=`ls -lrt $weeklyDir/* | awk '{print $NF}' | tail -n 1`
	dateFileLastBackup=`date -r $nameFileLastBackup  +"%s"` # date in (number of seconds since the epoch, when the seventies begun, UTC)
	diffInSec=`expr $dateFileCurrentBackup - $dateFileLastBackup`
	diffInDay=`expr $diffInSec / 86400 + $gap`

	# If diff with last backup is over this determined we copy it
	if (( $diffInDay >= $dayInterval ));then
		prepare_backup_to_store
        	# Move it in weekly directory
        	mv ${now}.tar.gz $weeklyDir/
	
		#Now check if we have to remove old backup
		nbBackup=`ls -lrt $weeklyDir/* | wc -l`
		if(($nbBackup > $nbWeeklyMax));then
			nameFileOlderBackup=`ls -lrt $weeklyDir/* | awk '{print $NF}' | head -n 1`	
			rm $nameFileOlderBackup
		fi
	fi
fi

###################
## Handle monthly backup
if [ ! "$(ls -A $monthlyDir)" ]; then # directory empty
        if [[ -f $weeklyDir/${now}.tar.gz ]];then #Current backup already managed, we just copy it
		cp $weeklyDir/${now}.tar.gz $monthlyDir/
	else # #Current backup not yet managed, we handle it
		prepare_backup_to_store
	        # Move it in montly directory
	        mv ${now}.tar.gz $monthlyDir/
	fi
else # directory not empty
	nameFileLastBackup=`ls -lrt $monthlyDir/* | awk '{print $NF}' | tail -n 1`
	dateFileLastBackup=`date -r $nameFileLastBackup  +"%s"` # date in (number of seconds since the epoch, when the seventies begun, UTC)
	diffInSec=`expr $dateFileCurrentBackup - $dateFileLastBackup + $gap`
        diffInMonth=`expr $diffInSec / 2629746`
	
	# If diff with last backup is over this	determined we copy it  	  
        if (( $diffInMonth >= $monthInterval ));then
		if [[ -f $weeklyDir/${now}.tar.gz ]];then #Current backup already managed, we just copy it
                	cp $weeklyDir/${now}.tar.gz $monthlyDir/
		else # #Current backup not yet managed, we handle it
			prepare_backup_to_store
                	# Move it in monthly directory
                	mv ${now}.tar.gz $monthlyDir/
        	fi

		#Now check if we have to remove old backup
		nbBackup=`ls -lrt $monthlyDir/* | wc -l`
	        if(($nbBackup > $nbMonthlyMax));then
                        nameFileOlderBackup=`ls -lrt $monthlyDir/* | awk '{print $NF}' | head -n 1`
                        rm $nameFileOlderBackup
                fi
	fi
fi


###################
## Handle yearly backup
if [ ! "$(ls -A $yearlyDir)" ]; then # directory empty
        if [[ -f $weeklyDir/${now}.tar.gz ]];then #Current backup already managed, we just copy it
                cp $weeklyDir/${now}.tar.gz $yearlyDir/
        else # #Current backup not yet managed, we handle it
		prepare_backup_to_store
                # Move it in montly directory
                mv ${now}.tar.gz $yearlyDir/
        fi
else # directory not empty
        nameFileLastBackup=`ls -lrt $yearlyDir/* | awk '{print $NF}' | tail -n 1`
        dateFileLastBackup=`date -r $nameFileLastBackup  +"%s"` # date in (number of seconds since the epoch, when the seventies begun, UTC)
        diffInSec=`expr $dateFileCurrentBackup - $dateFileLastBackup + $gap`
        diffInYear=`expr $diffInSec / 31556952`
	
        # If diff with last backup is over this determined we copy it
        if (( $diffInYear >= $yearInterval ));then
                if [[ -f $weeklyDir/${now}.tar.gz ]];then #Current backup already managed, we just copy it
                        cp $weeklyDir/${now}.tar.gz $yearlyDir/
                else # #Current backup not yet managed, we handle it
			prepare_backup_to_store
                        # Move it in monthly directory
                        mv ${now}.tar.gz $yearlyDir/
        	fi
	fi
fi

