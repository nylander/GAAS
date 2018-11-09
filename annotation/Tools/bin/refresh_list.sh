#!/bin/bash

#While looking at all script within the repo we skip all Deprecated folder and what it contains
for i in $(find ../.. -not \( -path */Deprecated -prune \) -name '*.pl' -o -name '*.sh' -o -name '*.py' -o -name '*.r' -o -name '*.R' -o -name '*.rb');do

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
