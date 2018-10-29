# CIGAR format overview

The CIGAR format is quite diverce and it is sometimes hard to understand it.
I found a really nice description of the history of the CIGAR format in the [LASTZ manual](http://www.bx.psu.edu/~rsharris/lastz/newer/README.lastz-1.02.40.html#ex_cigar)

I found other nice resources :
>[http://lh3.github.io/2018/03/27/the-history-the-cigar-x-operator-and-the-md-tag](http://lh3.github.io/2018/03/27/the-history-the-cigar-x-operator-and-the-md-tag)
>[https://github.com/vsbuffalo/devnotes/wiki/The-MD-Tag-in-BAM-Files](https://github.com/vsbuffalo/devnotes/wiki/The-MD-Tag-in-BAM-Files)
>[EnsemblDocs Wiki](http://htmlpreview.github.io/?https://github.com/NBISweden/GAAS/blob/master/annotation/CheatSheet/snapshots/ensembl_cigar.html)

**What we have to retain is:**

CIGAR is an acronym for Concise Idiosyncratic Gapped Alignment Report

## Original Exonerate CIGAR

CIGAR has been originally defined by the Exonerate alignment program.
 **Exonerate CIGAR format does not include nucleotides**. I didn't find the original description of the format from the exonerate web-pages (it has been removed), but we still can find it on other old ressources: </br>
from 2004 FlyBase here: [http://rice.bio.indiana.edu:7082/annot/gff3.html](http://rice.bio.indiana.edu:7082/annot/gff3.html) </br>
from 2010 WormBase here: [http://wiki.wormbase.org/index.php/GFF3specProposal](http://wiki.wormbase.org/index.php/GFF3specProposal) </br>

Here is what could have been the original Exonerate definition:

Operator | Description
-- | --
M    |    match
I     |   insert a gap into the reference sequence
D    |   insert a gap into the target (delete from reference)
F     |   frameshift forward in the reference sequence
R     |   frameshift reverse in the reference sequence

Each run is encoded by the letter code, whitespace, and the length; multiple runs are separated by whitespace. 
So it could looks like that:
>M 24 I 3 M 7 D 2 M 19

## Ensembl CIGAR
Then Ensembl has created the **ensembl cigar format** (didn't fdind any information about it until now),

## Samtools original CIGAR
**SAMtools created the extended cigar string**.
Here is a priori the original specification of the Samtools CIGAR format:

Operator | Description
-- | --
D | Deletion; the nucleotide is present in the reference but not in the read
H | Hard Clipping; the clipped nucleotides are not present in the read.
I | Insertion; the nucleotide is present in the read  but not in the reference.
M | Match; can be either an alignment match or mismatch. The nucleotide is present in the reference.
N | Skipped region; a region of nucleotides is not present in the read
P | Padding; padded area in the read and not in the reference
S | Soft Clipping;  the clipped nucleotides are present in the read

So it could looks like that:
>24M3I7M2D19M
In this variant, whitespace is removed and the order of the letter code and length are reversed (length appears before letter code).

## Samtools extended CIGAR

The “M” didn't allow to differentiate between the matches and mismatches. So, “X” (substitution/Mismatch) and “=”  (match) have been added to the specification in response to the request of several important users (c.f. [The history the MD tag and the CIGAR X operator](http://lh3.github.io/2018/03/27/the-history-the-cigar-x-operator-and-the-md-tag)). 

Here is the current specification of the Samtools CIGAR format:

Operator | Description
-- | --
D | Deletion; the nucleotide is present in the reference but not in the read
H | Hard Clipping; the clipped nucleotides are not present in the read.
I | Insertion; the nucleotide is present in the read  but not in the reference.
M | Match; can be either an alignment match or mismatch. The nucleotide is present in the reference.
N | Skipped region; a region of nucleotides is not present in the read
P | Padding; padded area in the read and not in the reference
S | Soft Clipping;  the clipped nucleotides are present in the read
X | Read Mismatch; the nucleotide is present in the reference
= | Read Match; the nucleotide is present in the reference

So it could looks like that:
>16=X7=3I7=2DX18=

/!\ In some variants the length is omitted if it is 1 (I don't know yet which version)

# New Exonerate CIGAR

if you go to the [exonerate manual web page](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate-manual) we can see that the specification has evolved a lot.
Here are their last description of the format. Some definition are now really clear and I didn't find any further definition (e.g 5 and 3).

Operator | Description
-- | --
M | Match
C | Codon
G | Gap
N | Non-equivalenced region
5 | 5' splice site
3 | 3' splice site
I | Intron
S | Split codon
F | Frameshift 

# All CIGAR operators gathered in one table

To gather all the operator information in one place,  I did a union of the different operators of the different formats and end-up with this last table:

Operator | Description
-- | --
M | Match ; can be either an alignment match or mismatch. The nucleotide is present in the reference.
C | Codon
G | Gap
N | Non-equivalenced region
5 | 5' splice site
3 | 3' splice site
I | Intron / the nucleotide is present in the read  but not in the reference. / insert a gap into the reference sequence
S | Split codon / Soft Clipping;  the clipped nucleotides are present in the read
H | Hard Clipping; the clipped nucleotides are not present in the read
F | Frameshift / frameshift forward in the reference sequence
D | Deletion; the nucleotide is present in the reference but not in the read / insert a gap into the target (delete from reference)
P | Padding; padded area in the read and not in the reference
X | Read Mismatch; the nucleotide is present in the reference
= | Read Match; the nucleotide is present in the reference
R     |   frameshift reverse in the reference sequence
