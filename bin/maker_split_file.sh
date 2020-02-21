#!/bin/bash

file=$1
script="~/git/code/Tools/Maker/split_gff_by_source.pl"
basename=$(basename $file .gff)

echo "Splitting file $file"
mkdir -p annotation/$basename

perl $script --input $file -d annotation/$basename

echo "Converting Maker file to GTF"
/sw/bioinfo/cufflinks-2.1.1/gffread -o annotation/$basename/maker.gtf -T -F annotation/$basename/maker.gff

echo "All done!"
