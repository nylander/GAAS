#!/bin/sh -e

# References:
#   https://support.10xgenomics.com/genome-exome/software/pipelines/latest/advanced/references
#   https://en.wikipedia.org/wiki/Bin_packing_problem

max_contigs=500
padding=500
out_width=70

order=

while getopts 'n:p:rs:w:' opt;
do
    case $opt in
        n) max_contigs=$OPTARG ;;
        p) padding=$OPTARG ;;
        r) order=r ;;
        w) out_width=$OPTARG ;;
        *)
            echo 'Error in command line parsing' >&2
            exit 1
    esac
done

shift "$(( OPTIND - 1 ))"

fasta_in="$1"
fasta_out="$(dirname "$fasta_in")/packed-$(basename "$fasta_in")"
fasta_map="$(dirname "$fasta_in")/$(basename "$fasta_in").map"

if [ ! -f "$fasta_in" ]; then
    printf 'Can not read from "%s"\n' "$fasta_in" >&2
    exit 1
fi

cat >&2 <<MESSAGE_END
Using the following parameters:

    max number of contigs:      (-n) $max_contigs
    padding:                    (-p) $padding
    output fasta width          (-w) $out_width
    output file:                     $fasta_out
MESSAGE_END
if [ "$order" = "r" ]; then
    echo '    sorting:                    (-r) decreasing' >&2
else
    echo '    sorting:                 (no -r) increasing' >&2
fi

padding_str=$( perl -e "print 'N' x $padding" )

samtools faidx "$fasta_in"

contigs=0

sort -k "2,2n$order" <"$fasta_in.fai" |
awk -f "$(dirname "$0")"/fa-firstfit.awk \
    -v padding="$padding" \
    -v maxcontigs="$max_contigs" |
sort -k1,1 -k3,3n |
tee "$fasta_map" |
while read contig seq len; do
    if [ "$contig" != "$prevcontig" ]; then
        prevcontig=$contig
        pad=0

        contigs=$(( contigs + 1 ))
        if [ "$contigs" -gt "$max_contigs" ]; then
            echo 'Error: Maximum number of contigs exceeded!' >&2
            exit 1
        fi
    else
        pad=1
    fi

    if [ "$pad" -eq 1 ]; then
        printf '%s\n' "$padding_str"
    else
        printf '>%s\n' "$contig"
    fi

    samtools faidx "$fasta_in" "$seq" | sed 1d
done |
awk -f "$(dirname "$0")"/fa-reformat.awk -v maxlen="$out_width" >"$fasta_out"
