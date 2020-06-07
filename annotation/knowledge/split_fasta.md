# Split a FASTA file
## Review of the main way to split a FASTA file


tool | language | One sequence per file | Can select chunck nb | Can select nb seq by chunck | Can select output file size | Overlap possible | Subsample possible | Example | Comment
-- | -- | -- | -- | -- | -- | -- | -- | -- | -- |
awk | awk | yes | probably | yes | no | probably | probably | [example](awk) | Not easy for novice 
split | bash | yes | no | yes | yes | no | no | [example](split) | Fasta must be single line fasta (one header + one single sequence line)
[AGAT](https://github.com/NBISweden/AGAT) | Perl | yes | yes | yes | no | yes | yes | [example](agat) | /
[PyFasta](https://pypi.org/project/pyfasta/#command-line-interface) | Python |  |  |  |  |  |  | [example](pyfasta) | 
[pyfaidx](https://github.com/mdshw5/pyfaidx) | Python |  |  |  |  |  |  | [example](pyfaidx) |

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

## AGAT

## PyFasta

split a fasta file into 6 new files of relatively even size:

`pyfasta split -n 6 original.fasta`

## pyfaidx

`faidx --split-files original.fasta`

## Interesting thread

https://www.biostars.org/p/229441/
