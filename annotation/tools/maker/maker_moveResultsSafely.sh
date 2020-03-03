#!/bin/bash

############################################################################
# JD 2014/04
# What it does ? clean move of maker results
# use:  script Directory 
############################################################################

# Arguments and Paths
logFile=""
if (( $# !=1 )); then
        echo -e "You have to specify a new directory ..."
	exit
else
	dirRes=$1
	if [ ! -d $dirRes ];then
		mkdir $dirRes
	else
		echo -e "The directory already exists !\nDo you Want Overwrite this directory ?\n	(0) yes\n	(1) no\n	(2) Take only things absent in the target directory"
		read
            	case $REPLY in
                        0) echo "Let s go to delete the directory"
	                	rm -r $dirRes        
				mkdir $dirRes
			;;
                        1) echo "Ok we keep it ! See you then."
                        ;;
			2) echo "I will move only directories/Files that are absent in the target directory"
			moveifabsent="yes"
			;;
                esac

	fi
fi


# Directory move
if [ ! -d genome.maker.output ];then
	echo "/!\\ <genome.maker.output> folder not found !"
else
	if [[ moveifabsent == "yes" ]];then
		if [ ! -d $dirRes/genome.maker.output ];then
			mv genome.maker.output $dirRes/
		fi
	else
		mv genome.maker.output $dirRes/
	fi
fi

if [ ! -d annotation ];then
        echo "/!\\ <annotation> folder not found !"
	echo "Suggestion: Did you think to launch the script to perform this directory ?"
else
    	if [[ moveifabsent == "yes" ]];then
                if [ ! -d $dirRes/annotation ];then
                        mv annotation $dirRes/
                fi
        else
                mv annotation $dirRes/
        fi
fi



#Files Move
if [ ! -f annotations.gff ];then
        echo "/!\\ <annotations.gff> file not found !"
	echo "Suggestion: Did you think to launch the script to perform this file ?"
else
	if [[ moveifabsent == "yes" ]];then
		if [ ! -f $dirRes/annotations.gff ];then
			mv annotations.gff $dirRes/
		fi
	else
	    	mv annotations.gff $dirRes/
	fi
fi

if [ ! -f annotations.proteins.fa ];then
        echo "/!\\ <annotations.proteins.fa> file not found !"
else
        if [[ moveifabsent == "yes" ]];then
                if [ ! -f $dirRes/annotations.proteins.fa ];then
                        mv annotations.proteins.fa $dirRes/
                fi
        else
    		mv annotations.proteins.fa $dirRes/
	fi
fi

if [ ! -f maker_opts.ctl ];then
        echo "/!\\ <maker_opts.ctl> file not found !"
else
	if [[ moveifabsent == "yes" ]];then
                if [ ! -f $dirRes/maker_opts.ctl ];then
                        cp maker_opts.ctl $dirRes/
                fi
        else
                cp maker_opts.ctl $dirRes/
        fi
    cp maker_opts.ctl $dirRes/
fi
