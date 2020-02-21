#!/bin/bash

        #remove scripts
        for j in *;do
                name=$(basename $j)
                if [ ! $name == "gaas_refresh_list.sh"  ] ; then
                        unlink $j
                fi
        done

#While looking at all script within the repo we skip all Deprecated folder and what it contains
for i in $(find .. -not \( -path */Deprecated -prune \) -name '*.pl' -o -name '*.sh' -o -name '*.py' -o -name '*.r' -o -name '*.R' -o -name '*.rb');do

        name=$(basename $i)

        #Populate scripts
        #if script does not exist, create a link
        if [[ ! -f $name ]];then
                ln $i
        fi

done
