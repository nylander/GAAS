# Utility for preparing an assembly for polishing using 10X Chromium data

Authors:
Original code by Andreas Kähäri, modified by Mahesh Binzer-Panchal. 

## Description

Polishing with 10X Chromium data in theory feels like an easy task. However, in order
to use the 10X tool `longranger align`, the reference (assembly) must be formatted to
have certain properties. These are:

* Diploid genome — The phasing algorithm assumes 2 haplotypes. (Not relevant for the 
`align` submodule)
* 500 FASTA entries or fewer — if your assembly has more than 500 FASTA entries, 
concatenate smaller contigs together with 500 N's separating each original contig, 
until there are fewer than 500 FASTA entries total. Note: when creating a concatenated 
reference contig, you must create a `primary_contigs.txt` file, and omit the 
concatenated contig from `primary_contigs.txt`. See below for details.
* All contigs must be no more than 2^29-1 bp, or 528Mb, in length; this is a limitation 
of the BAM index file format.
* All contigs must have no colons or spaces in their names.

This is described on the 10X genomics website (https://support.10xgenomics.com/genome-exome/software/pipelines/latest/advanced/references)
as how to build a reference. The resulting alignment then includes trimmed 10X
genomics reads aligned in a GEM barcode aware manner. 

A note on 10X raw data.
* Adapter read through is present in raw data.
* The first 16bp of Read 1 is the 10X barcode. If the Illumina run was restarted for
any reason, the first base of the GEM barcode may be missing causing issues, but the 
sequencing provider should know to do a complete restart.
* The filename must be formatted in the same way `longranger mkfastq` would format it:
`Sample1_S1_L001_I1_001.fastq.gz` (Index file - May not be necessary), `Sample1_S1_L001_R1_001.fastq.gz` (Read 1), 
`Sample1_S1_L001_R2_001.fastq.gz` (Read 2).

We have developed a script that merges the smaller contigs into larger psuedoscaffolds,
while trying to keep the largest unconcatenated.

## How to use

Given a diploid assembly, merge the two assemblies. It is important that the assembly
does not contain N's, or this will complicate unpacking later:

```bash
set -ue
PRIMARY_ASSEMBLY=primary_assembly.fasta
ALTERNATE_ASSEMBLY=alt_assembly.fasta
MERGED_ASSEMBLY=merged_assembly.fasta
cat "$PRIMARY_ASSEMBLY" "$ALTERNATE_ASSEMBLY" > "$MERGED_ASSEMBLY"
# Check N count
grep -v ">" "$MERGED_ASSEMBLY" | tr -dc "Nn" | wc -m
```

Then pack the merged assembly, make a `primary_contigs.txt` file, and a reference:

```bash
binpack/prep-contigs.sh "MERGED_ASSEMBLY" &> "${MERGED_ASSEMBLY}.pack_log"
# Use grep to make a primary_contigs.txt file - Here I use the keyword arrow since concatenated contigs do not have this
grep arrow "${MERGED_ASSEMBLY}.pack_log" > primary_contigs.txt
# Make a reference
longranger mkref "packed-${MERGED_ASSEMBLY}"
REFERENCE="refdata-packed-${MERGED_ASSEMBLY/.fasta}"
```

The 10X data can then be aligned using `longranger align`:

```bash
longranger align --id='some_name' --reference="$REFERENCE" --fastqs='/path/to/fastqs/directory'
```

The output can be found in `some_name/outs`. This can then be polished using Pilon (Pilon 
estimates memory needs of 1GB per mb of sequence - The java flag `-d64` runs java in 64bit
mode allowing for higher specifications of heap memory `Xmx`):

```bash
BAMFILE=some_name/outs/possorted_bam.bam
GENOME="$REFERENCE/fasta/genome.fa"
PREFIX='output_prefix'
java -d64 -Xmx2t -jar pilon-1.23.jar --threads "$CPUS" \
   --genome "$GENOME" --frags "$BAMFILE" \
   --outdir "${PREFIX}-polished" \
   --output "${PREFIX}" --tracks --changes
```

You should now have a polished genome, IGV viewable tracks, and changes made to the packed assembly.

During packing, a list was made of which contigs were put into which packed contigs. We now restore
the old contig names while unpacking the sequence (splitting on N's), and appending `|pilon` to note 
the contigs were polished using Pilon:

```bash
# The map file was created in the packing process
MAPFILE="$MERGED_ASSEMBLY.map"
POLISHED_GENOME="${PREFIX}-polished/${PREFIX}.fasta"
awk -v repl="$MAPFILE" -f binpack/unpack_contigs.awk "$POLISHED_GENOME" > "${PREFIX}_unpacked_polished.fasta"
```

## Other notes

### Use of Racon as an alternative to Pilon

Racon is another tool that can be used for polishing, that is supposed to be better and less resource 
intensive for larger assemblies. However, no list of changes can be created as output directly from
the tool, and furthermore more segments of N's appear in the polished assembly (either created by changing 
known bases or correcting inside a padding sequence - I don't know). This means unpacking would
require a dynamic programming algorithm to correctly unpack the corrected contigs and correctly assign
their original names. 

This is how one can use racon to polish the genome:

```bash
# Samtools, parallel, and racon are added to the path.
export PATH="$PATH:tools/samtools-1.9"
export PATH="$PATH:tools/parallel-20181022/src"
export PATH="$PATH:tools/racon/build/bin"

CPUS="${SLURM_NPROCS:-10}"

## Rename reads for compatibility with racon

TMP_READDIR="${SNIC_TMP:-/tmp}/reads"
mkdir -p "$TMP_READDIR"

BAMFILE=some_name/outs/possorted_bam.bam
TMP_SAMFILE=$(mktemp -p "$TMP_READDIR" --suffix ".sam")
TMP_ALLREADS=$(mktemp -p "$TMP_READDIR" --suffix ".fastq")
samtools view -@ "$CPUS" -h "$BAMFILE" | parallel --pipe awk -f mod_SAM.awk | tee "$TMP_SAMFILE" | parallel --pipe awk -f SAM2FASTQ.awk > "$TMP_ALLREADS"

## Run racon with updated fastqs and samfile.

GENOME="$REFERENCE/fasta/genome.fa"
PREFIX='output_prefix'

racon -u -t "$CPUS" "$TMP_ALLREADS" "$TMP_SAMFILE" "$GENOME" > "${PREFIX}_polished_racon.fasta"

rm -r "$TMP_READDIR"

```

where `mod_SAM.awk` is:

```awk
BEGIN {
        OFS="\t"
} 

{
        if(!/^@/) {
                if(and($2,0x40)){ 
                        $1=$1":1"
                } else if (and($2,0x80)) {
                        $1=$1":2"
                } 
        }
        print $0 
}
```

and `SAM2FASTQ.awk` is (this could probably be replaced with `samtools fastq`):

```awk
!/^@/ && !and($2,0x100) { 
        print "@"$1"\n"$10"\n+\n"$11 
}
```
