#!/bin/bash

NAME=$1

if [ -z "$NAME" ]
then
	echo "Must provide a name!"
else
	fathom -categorize 1000 genome.ann genome.dna
	fathom -export 1000 -plus uni.ann uni.dna
	forge export.ann export.dna
	hmm-assembler.pl $NAME . > $NAME.hmm
fi
