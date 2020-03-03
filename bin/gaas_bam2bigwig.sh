#!/bin/bash
# Convert BAM to BigWig
# Requires genome.list (tab-delimited list of chromosomes and lengths)

file=$1

# load bamToBed
module load bedtools/2.17.0
# load bedGraphToBigWig 
module load ucsc-tools

command -v bamToBed >/dev/null 2>&1 || { echo >&2 "I require bamToBed but it's not in PATH.  Aborting."; exit 1; }
command -v bedGraphToBigWig >/dev/null 2>&1 || { echo >&2 "I require bedGraphToBigWig but it's not in PATH.  Aborting."; exit 1; }

echo "Converting file to BED"
bamToBed -i $file > $(basename $file .bam).bed
echo "Creating coverage track"
genomeCoverageBed -i $(basename $file .bam).bed -bg -g genome.list > $(basename $file .bam).cov
echo "Writing BigWig file"
bedGraphToBigWig $(basename $file .bam).cov genome.list $(basename $file .bam).bw
