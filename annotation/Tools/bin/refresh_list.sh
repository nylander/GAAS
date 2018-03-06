#!/bin/bash

for i in $(find ~/git/NBIS/GAAS/annotation/ -name '*.pl' -o -name '*.sh' -o -name '*.py' -o -name '*.r' -o -name '*.R');do

	name=$(basename $i)
	
	#Populate scripts
	#if script does not exist, create a link
	if [[ ! -f $name ]];then
		ln -s $i 
	fi 
	
	#remove script no existing anymore
	# test if file exists (test actual file, not symbolic link)
	for j in *;do
		if [ ! -e "$j" ] ; then
			echo "$j link broken"
    		 	# code if the symlink is broken
			unlink $j
		fi
	done

	
done
