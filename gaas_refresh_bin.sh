#!/bin/bash

cd bin

if [[ ! ${PWD##*/} == "bin" ]];then
	echo "The script must be run in GAAS/bin folder"
	echo "currently here `pwd`"
	exit
fi


#remove scripts
for j in *;do
        name=$(basename $j)
        unlink $name
done

#While looking at all script within the repo we skip all Deprecated folder and what it contains
for i in $(find ../ -not \( -path */Deprecated -prune \) -not \( -path */bin -prune \)  -not \( -path */blib -prune \) -name '*.pl' -o -name '*.sh' -o -name '*.py' -o -name '*.r' -o -name '*.R' -o -name '*.rb');do

        name=$(basename $i)

	# skip gaas_refresh_list.sh because must not be in the bin to avoid to be distributed
	if [ $name == "gaas_refresh_bin.sh"  ] ; then
		continue
	fi

        #Populate scripts using hard link
        if [[ ! -f ${name} ]];then
                ln $i ${name}
        fi

done
