# CIGAR format overview

The CIGAR format is quite diverce and it is sometimes hard to understand it. Here is an overview of the format and its history.

**Important ressources that has been used to write this overview:**

>[http://lh3.github.io/2018/03/27/the-history-the-cigar-x-operator-and-the-md-tag](http://lh3.github.io/2018/03/27/the-history-the-cigar-x-operator-and-the-md-tag)  
>[https://github.com/vsbuffalo/devnotes/wiki/The-MD-Tag-in-BAM-Files](https://github.com/vsbuffalo/devnotes/wiki/The-MD-Tag-in-BAM-Files)  
>[EnsemblDocs Wiki](http://htmlpreview.github.io/?https://github.com/NBISweden/GAAS/blob/master/annotation/CheatSheet/snapshots/ensembl_cigar.html)  
[LASTZ manual](http://www.bx.psu.edu/~rsharris/lastz/newer/README.lastz-1.02.40.html#ex_cigar)  

**Forewords:**

CIGAR is an acronym for **C**oncise **I**diosyncratic **G**apped **A**lignment **R**eport and has been originally defined by the **Exonerate** alignment program. It is designed to contain the minimal information necessary for the reconstruction of an alignment. One alignment is described per line, to allow easy manipulation with UNIX tools.**Exonerate CIGAR format does not include nucleotides**.  
It exists other related format: 
 **Sugar** - Simple Ungapped Alignment Report  
 **Vulgar** - Verbose Ugly Labelled Gapped Alignment Report  

## Original Exonerate CIGAR (~2003)

Cigar format looks like this:

`cigar: hs989235.cds 5 468 + hsnfg9.embl 25689 27450 + 1916 M 13 I 1 M 35 I 1 M 4 I 1 M 13 D 1 M 4 I 1 M 115 D 404 M 37 D 1 M 164 I 1 M 12 D 898 M 16 I 1 M 12 I 1 M 21 D 1 M 10`

The fields are as follows:

>1. query identifier  
>2. query start position  
>3. query stop position  
>4. query strand  
>5. target identifier  
>6. target start position  
>7. target stop position  
>8. target strand  
>9. score  
>10. **CIGAR string**

The 10th field is often called **CIGAR string** and are defined in pairs. Each pair is also called a run. The cigar string describes the edit path throught the alignment. These contain a M,I,D or N corresponding to a Match, Insert, Delete or iNtron, followed by the length.

Operator | Description
-- | --
M    |    **M**atch
I     |   **I**nsert
D    |   **D**elete
N     |   i**N**tron
 
/!\ I havn't find any ressource using/talking about the N (iNtron) operator. It seems to be used by D suystemtically. Any information about it is very welcome.

To resume:  
Each pair/run is encoded by the letter code/operator, whitespace, and the length; multiple pairs/runs are separated by whitespace.

## Ensembl CIGAR
Then Ensembl has created the **ensembl cigar format**.

In the Ensembl CIGAR format the numbers and letters are switched, and there are no gaps in the string. So the above example in Ensembl would appear in a feature table in three rows with these CIGAR strings:

`>13M1I35M1I4M1I13M1D4M1I115M  
>37M1D164M1I12M  
>16M1I12M1I21M1D10M`

In the Ensembl CIGAR format the numbers and letters are switched, and there are no gaps in the string. So the above example in Ensembl would appear in a feature table in three rows with these CIGAR strings:

## Updated Exonerate CIGAR - Gap attribute in GFF3(~2004-2005)
 
 [ref](https://doi.org/10.1186/1471-2105-6-31)  
I didn't find the original description of the format from the exonerate web-pages (it has been removed), but we still can find it on other old ressources: </br>
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
