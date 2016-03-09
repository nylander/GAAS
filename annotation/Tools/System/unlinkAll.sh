#!/bin/bash

# JD 2014
#
#The script allow to unlink all links in a directory
####



# Arguments and Paths
############################################################################

if (( $# != 1 )) ; then
        echo -e "The script needs 1 parameter: \n(1)Directory where all the link will be unlinked"
        exit
fi

IFS=$'\n'
for i in $(ls -l $1);do
        if [[ $i =~ ^[l] ]];then
                paired=$(echo $i | awk '{print $9" "$11}')
                to=$(echo $paired | cut -d' ' -f1)
                unlink $1/$to
        fi
done

