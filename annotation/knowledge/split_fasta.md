# Split a FASTA file
## Review of the main way to split a FASTA file


tool | language | One sequence per file | Can select chunck nb | Can select nb seq by chunck | Can select output file size | Overlap possible | Can cut sequence | Subsample possible | Example | Comment
-- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- |
awk | awk | yes | no | yes | no | no | no | no | [example](awk) | Not easy for novice 
split | bash | yes | no | yes | yes | no | no | no | [example](split) | Fasta must be single line fasta (one header + one single sequence line)
bash | bash | yes | no | no | no | no | no | no | [example](bash) |  Individual files will have the name of the corresponding sequence, without leading >
gaas_fasta_splitter.pl from [GAAS](https://github.com/NBISweden/GAAS) | Perl | yes | yes | yes | no | yes | yes | yes | [example](agat) | /
[PyFasta](https://pypi.org/project/pyfasta/#command-line-interface) | Python |  |  |  |  |  |  |  | [example](pyfasta) | 
[pyfaidx](https://github.com/mdshw5/pyfaidx) | Python |  |  |  |  |  |  |  | [example](pyfaidx) |
[GenomeTools](https://github.com/genometools/genometools) |  Mostly C | yes | yes | no | yes | no | no | no | [example](GenomeTools) | 
[seqretsplit](http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/seqretsplit.html) from [EMBOSS](http://emboss.sourceforge.net/what/) |  C | yes | no | no | no | no | no | no | [example](EMBOSS) |
bp_seqretsplit.pl from [Bioperl](https://github.com/bioperl/bioperl-live) |  perl | yes | no | no | no | no | no | no | [example](bp_seqretsplit) |


# Example

## Awk

size = chunk size
pre = output file prefix
pad = padding width (the width of the numeric suffix).

```
awk -v size=1000 -v pre=prefix -v pad=5 '
   /^>/ { n++; if (n % size == 1) { close(fname); fname = sprintf("%s.%0" pad "d", pre, n) } }
   { print >> fname }
' input.fasta
```

## Split

`split -l 2000 input.fasta`

## Bash

```
while read line
do
    if [[ ${line:0:1} == '>' ]]
    then
        outfile=${line#>}.fa
        echo $line > $outfile
    else
        echo $line >> $outfile
    fi
done < myseq.fa
```

## GAAS

`gaas_fasta_splitter.pl  --nb_chunks 10 --nb_seq_by_chunk 100 --overlap 50 --size_seq 1000000`

## PyFasta

split a fasta file into 6 new files of relatively even size:

`pyfasta split -n 6 original.fasta`

## pyfaidx

`faidx --split-files original.fasta`

## GenomeTools

`gt splitfasta -splitdesc multifastafile.fa`

## EMBOSS

`seqretsplit input.fa`

## bp_seqretsplit

`bp_seqretsplit file1 file2`

Similar to:  
```
#!/usr/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;
my $in = Bio::SeqIO->new(-format => 'fasta',
                         -fh   => \*ARGV);
while( my $s = $in->next_seq ) {
    my ($id) = ($s->id =~ /^(?:\w+)\|(\S+)\|/);
    Bio::SeqIO->new(-format => 'fasta',
                    -file   => ">".$id.".fasta")->write_seq($s);
}
```

## Reference

https://www.biostars.org/p/229441/
