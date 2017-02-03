Let's develop a tool to detect motochondrial contigs within assemblies

# MitoSearch

We use min hash sketches of the Mitochonrial and Plasmid RefSeq sequences stored
in a Sequence Bloom Tree (SBT) to filter and probablistically match if contigs contain
mitochondrial or plasmid sequences. 

If a contig does not have a close relation to a min hash signature in the tree, then
the sequence is definitely not a mitochondrial or plasmid sequence in the database.
If a match is found, it is highly likely that the sequence is mitochondrial or plasmid, 
but not definite.

We use the SBT and min hash features implemented in SourMash (https://github.com/dib-lab/sourmash.git)
to do this search. GNU parallel can be used to speed up searching. SeqTK is used for convenience.

## Installation

See INSTALL.md

## Usage

### Quick Start

A Sequence Bloom Tree of the Mitochondrial and Plasmid RefSeq sequences is provided
in the SBT folder. 
All one needs to do is create min hash signatures of their contigs, and then search
the SBT.

```
# Separate all contigs from your draft assembly into their own files.
seqtk seq -A draft_assembly.fasta.gz | awk 'BEGIN { RS=">" } $0 != "" { print ">"$0 > $1".fasta" }'
# Compute sourmash min hash signatures of all those contigs.
sourmash compute -k21,31,51 -f *.fasta
# Search each of those signatures against the SBT to see if a match is found
for CONTIG_SIGNATURE in *.fasta.sig ; do
	echo $CONTIG_SIGNATURE
	sourmash sbt_search organelle $CONTIG_SIGNATURE
done > mitosearch.txt
# Find the contigs that have a match to the organelle database
grep -B2 ".fasta$" mitosearch.txt | less -S
# Find the organelles that have been matched
for MATCH in $(grep ".fasta$" mitosearch.txt | cut -f2 -d" " | sort -u ); do 
	grep ">" 0R_Reference/$MATCH
done
```

### Building the organelle SBT database

Here are the instructions on how to build the Organelle database (or database of your choice)

```
mkdir 0R_Reference
cd 0R_Reference
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.1.genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.2.1.genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/plastid.1.1.genomic.fna.gz
# Save individual organelle sequences in files named by their accession header (this could be sped up using parallel too)
zcat *.fna.gz | awk 'BEGIN { RS = ">" } $0 != "" { filename=$1; gsub(/\|/,"_",filename); print ">"$0 > filename".fasta" }'
# Compute the signatures of all the organelle sequences using k-mer sizes of 21, 31 (def.), and 51.
parallel -j <cores> sourmash compute -k21,31,51 -f -o {}.sig {} ::: *.fasta
# Make the SBT 
sourmash sbt_index organelle *.sig
```

### Searching the organelle SBT database

This is covered in the Quickstart but here are the same instructions using GNU parallel to speed things up.

```
zcat draft_assembly.fasta.gz | awk 'BEGIN { RS = ">" } $0 != "" { filename=$1; gsub(/\|/,"_",filename); print ">"$0 > filename".fasta" }'
parallel -j <cores> sourmash compute -k21,31,51 -f -o {}.sig {} ::: *.fasta
parallel -j <cores> sourmash sbt_search organelle {} ::: *.fasta.sig > mitosearch.txt
grep -B2 ".fasta$" mitosearch.txt | less -S
grep ".fasta$" mitosearch.txt | cut -f2 -d" " | sort -u | parallel -j <cores> grep ">" 0R_Reference/{}
```

### Blasting matches to organelle sequence

Although a match is shown to an organelle, blasting the contigs will provide additional confirmation of the contig origin.

```
# Get the contigs that have matches
grep -B2 ".fasta$" mitosearch.txt | grep ".fasta.sig" | sed -e 's/.sig//g' > mito_tigs.fofn
# Print the statistics of the contigs (Output format: contig, length, #A, #C, #G, #T, #2, #3, #4, #CpG, #tv, #ts, #CpG-ts)
parallel -a mito_tigs.fofn seqtk comp {} 
# Blast the individual contigs
parallel -a mito_tigs.fofn -j <cores> 'blastn -db $BLASTDB/nt -query {} -outfmt "6 qseqid sseqid stitle evalue" -max_hsps 1 -num_threads <threads> -evalue <evalue> | head -n1' | tee blast.log
```
On Milou use the following to load the necessary modules.
```
module load bioinfo-tools blast/2.5.0+ gnuparallel/20150522 seqtk/1.0-r68e
```
