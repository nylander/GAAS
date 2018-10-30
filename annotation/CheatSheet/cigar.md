# CIGAR format overview

The CIGAR format is quite diverce and it is sometimes hard to understand it. Here is an overview of the format and its history.

### Important resources that has been used to write this overview:

>[http://lh3.github.io/2018/03/27/the-history-the-cigar-x-operator-and-the-md-tag](http://lh3.github.io/2018/03/27/the-history-the-cigar-x-operator-and-the-md-tag)  
>[https://github.com/vsbuffalo/devnotes/wiki/The-MD-Tag-in-BAM-Files](https://github.com/vsbuffalo/devnotes/wiki/The-MD-Tag-in-BAM-Files)  
>[EnsemblDocs Wiki](http://htmlpreview.github.io/?https://github.com/NBISweden/GAAS/blob/master/annotation/CheatSheet/snapshots/ensembl_cigar.html)  
>[LASTZ manual](http://www.bx.psu.edu/~rsharris/lastz/newer/README.lastz-1.02.40.html#ex_cigar)  
>[http://rice.bio.indiana.edu:7082/annot/gff3.html](http://rice.bio.indiana.edu:7082/annot/gff3.html)  
>[http://wiki.wormbase.org/index.php/GFF3specProposal](http://wiki.wormbase.org/index.php/GFF3specProposal)  
>[exonerate manual web page](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate-manual)  

### Forewords:

CIGAR is an acronym for **C**oncise **I**diosyncratic **G**apped **A**lignment **R**eport and has been originally defined by the **Exonerate** alignment program. 
The CIGAR format is designed to contain the minimal information necessary for the reconstruction of an alignment. One alignment is described per line, to allow easy manipulation with UNIX tools.**Exonerate CIGAR format does not include nucleotides**.  

It exists other related format:    
  - **Sugar** - Simple Ungapped Alignment Report  
  - **Vulgar** - Verbose Ugly Labelled Gapped Alignment Report  
 
**What means idiosyncratic ?**  
  An idiosyncrasy is an unusual feature of the tool (though there are also other uses). It also means odd habit. The term is often used to express eccentricity or peculiarity.
Below is the idiodyncracies/conventions describefd in the man page of exonerate-1.0.0:


    CONVENTIONS
       A  number  of conventions (and idiosyncracies) are used within exonerate.  An understanding of them facili-
       tates interpretation of the output.

       Coordinates
              An in-between coordinate system is used, where the positions are counted between the symbols, rather
              than on the symbols.  This numbering scheme starts from zero.  This numbering is shown below for the
              sequence "ACGT":

               A C G T
              0 1 2 3 4

              Hence the subsequence "CG" would have start=1, end=3, and length=2.  This coordinate system is  used
              internally  in  exonerate,  and for all the output formats produced with the exception of the "human
              readable" alignment display and the GFF output where convention and standards dictate otherwise.

       Reverse Complements
              When an alignment is reported on the reverse complement of a sequence, the  coordinates  are  simply
              given  on  the  reverse complement copy of the sequence.  Hence positions on the sequences are never
              negative.  Generally, the forward strand is indicated by '+', the reverse  strand  by  '-',  and  an
              unknown or not-applicable strand (as in the case of a protein sequence) is indicated by '.'

       Alignment Scores
              Currently,  only  the raw alignment scores are displayed.  This score just is the sum of transistion
              scores used in the dynamic programming.  For example, in the case  of  a  Smith-Waterman  alignment,
              this will be the sum of the substitution matrix scores and the gap penalties.


## Original Exonerate CIGAR (before 2003)

### The related SUGAR format
We cannot talk about the CIGAR format without talking first about the **SUGAR** format. Sugar is Simple UnGapped Alignment Report, which displays ungapped alignments one-per-line. The sugar line starts with the string "sugar:" for easy extraction from the output, and is followed by the the following 9 fields in the order below:  

>1. query identifier  
>2. query start position  
>3. query stop position  
>4. query strand  
>5. target identifier  
>6. target start position  
>7. target stop position  
>8. target strand  
>9. score  

Here an example of SUGAR format:  
`sugar: hs989235.cds 5 468 + hsnfg9.embl 25689 27450 + 1916`

### The CIGAR format
The  format starts  with  the same 9 fields as sugar output (see above), and is followed by a series of <operation, length> pairs where operation is one of match, insert or delete, and the length describes the number of times this operation is repeated
Cigar format looks like this:

`cigar: hs989235.cds 5 468 + hsnfg9.embl 25689 27450 + 1916 M 13 I 1 M 35 I 1 M 4 I 1 M 13 D 1 M 4 I 1 M 115 D 404 M 37 D 1 M 164 I 1 M 12 D 898 M 16 I 1 M 12 I 1 M 21 D 1 M 10`

The 10th field is actually the **CIGAR string** and are defined in pairs. Each pair is also called a run. The cigar string describes the edit path throught the alignment. These contain a M,I or D corresponding to a Match, Insert or Delete, followed by the length.

Operator | Description
-- | --
M    |    **M**atch
I     |   **I**nsert
D    |   **D**elete


To resume:  
Each pair/run is encoded by the letter code/operator, whitespace, and the length; multiple pairs/runs are separated by whitespace.

## Ensembl CIGAR
Then Ensembl has created the **ensembl cigar format**.

In the Ensembl CIGAR format the numbers and letters are switched, and there are no gaps in the string. So the above example in Ensembl would appear in a feature table in three rows with these CIGAR strings:

>13M1I35M1I4M1I13M1D4M1I115M  
>37M1D164M1I12M  
>16M1I12M1I21M1D10M

/!\ What about N for i**N**tron ?
The ensembl page describe the CIGAR string from exonerate like that:  
`The cigar string describes the edit path throught the alignment. These contain a M,I,D or N corresponding to a Match, Insert, Delete or iNtron, followed by the length.`
I don't know where that N is coming from, I havn't find any ressource using/talking about the N (iNtron) operator. It seems to be used by D systemtically. Any information about it is very welcome.

## Updated Exonerate CIGAR - Gap attribute in GFF3(~2004)
 
The GFF3 format integrated a Gap attribute in the 9th column of the gff" files to describe alignements. Here is what has been described:



    Gap   The alignment of the feature to the target if the two are
          not colinear (e.g. contain gaps).  The alignment format is
	  taken from the CIGAR format described in the 
	  Exonerate documentation.
	  (http://cvsweb.sanger.ac.uk/cgi-bin/cvsweb.cgi/exonerate
           ?cvsroot=Ensembl).  See "THE GAP ATTRIBUTE" for a description
	   of this format.

Here is what could have been the original Exonerate definition:

Operator | Description
-- | --
M    |    match
I     |   insert a gap into the reference sequence
D    |   insert a gap into the target (delete from reference)
F     |   frameshift forward in the reference sequence
R     |   frameshift reverse in the reference sequence

As the original CIGAR string from Exonerate each run is encoded by the letter code, whitespace, and the length; multiple runs are separated by whitespace. 
So it could looks like that:
>M 24 I 3 M 7 D 2 M 19

/!\ For segments of length 1 the number can be omitted, so "8M1D6M" is equal to "8MD6M". 

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

In this variant, as for the description of the Ensembl CIGAR format, whitespace is removed and the order of the letter code and length are reversed (length appears before letter code) compared to the Exonerate CIGAR string. At the difference of the Ensembl CIGAR string, everything is encoded in only one line.
 So it could looks like that:
>24M3I7M2D19M

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

# The Exonerate VULGAR format

Here is the last specification of the format. Some definition are not really clear and I didn't find any further definition (e.g 5 and 3). 

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
