The gene-finding format / general feature format  
GFF / GFF1 / GFF2 / GFF2.5 / GFF3 / GTF / GTF2.2 formats
===========================

It's often hard to understand and differentiate all GFF/GTF formats/flavors. Here is an overview of the format and its history to help disentangle this complexity.


## Forewords:
⇨	When I use the term gff it includes all gff formats/flavors. (The first version of the format was not called gff1 but gff. But to make it easier I will always express the version of the format saying gff1 when it's the first version of it. So from now when I say gff it means all gff formats/flavors).  
⇨	I the same way, when I use the term gtf it includes all gtf formats/flavors.  
⇨	I have created the term **gxf** that means all the gff and gtf formats/flavors.

The GFF Protocol Specification was initially proposed by Richard Durbin and David Haussler.
The format has been originaly developed to help the gene prediction (or gene finding) world. Indeed the gene prediction methods are based on two main steps, **first** finding signals (starts, splice sites, stops, motifs, etc.) and regions (exons, introns, protein domains etc.); **secondly** combining these to give complete gene, RNA transcript or protein structures. These two steps were usually performed within the same program. In order to decoupled them they have created the format called GFF ('Gene-Finding Format') allowing the transfer of feature information from a tool to another one.

The GFF fomat has been developed to be easy to parse and process by a variety of programs in different languages (e.g Unix tools as grep and sort, perl, awk, etc). For these reasons, they decided that each feature is described on a single line, and line order is not relevant.

A GFF record is an extension of a basic (name,start,end) tuple (or "NSE") that can be used to identify a substring of a biological sequence. 

## GFF0(before 1997-11-13)

There is no clear information about how look the format at that time but it was close to the GFF1 format specification without the field "source" added the 1997-11-13.

## GFF1(1997-11-13):
Here the oldest complete description I found of the GFF1 format: https://web.archive.org/web/19980222142332/http://www.sanger.ac.uk:80/~rd/gff.html
[Here a snapshot of the olderst description of the format I found (2000)]((http://htmlpreview.github.io/?https://github.com/NBISweden/GAAS/blob/master/annotation/CheatSheet/snapshots/GFF_Spec.html).

I consider the format as GFF1 when they definitly defined the 9 fields of the format (1997-11-13 rd: added extra "source" field as discussed at Newton Institute meeting 971029). Before that the format was existing but was at the stage of version 0.

This GFF1 format contains 8 madatory fields and 9th one optional. Fields are:  

    <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [group]

Definition of these fields are:

    <seqname>
        The name of the sequence. Having an explicit sequence name allows a feature file to be prepared for a data set of multiple sequences. Normally the seqname will be the identifier of the sequence in an accompanying fasta format file. An alternative is that 'seqname' is the identifier for a sequence in a public database, such as an EMBL/Genbank/DDBJ accession number. Which is the case, and which file or database to use, should be explained in accompanying information.
    <source>
        The source of this feature. This field will normally be used to indicate the program making the prediction, or if it comes from public database annotation, or is experimentally verified, etc.
    <feature>
        The feature type name. We hope to suggest a standard set of features, to facilitate import/export, comparison etc.. Of course, people are free to define new ones as needed. For example, Genie splice detectors account for a region of DNA, and multiple detectors may be available for the same site, as shown above.
    <start>, <end>
        Integers. <start> must be less than or equal to <end>. Sequence numbering starts at 1, so these numbers should be between 1 and the length of the relevant sequence, inclusive.
    <score>
        A floating point value. When there is no score (i.e. for a sensor that just records the possible presence of a signal, as "splice5" above) you must give something, by convention 0.
    <strand>
        One of '+', '-' or '.'. '.' should be used when strand is not relevant, e.g. for dinucleotide repeats.
    <frame>
        One of '0', '1', '2' or '.'. '0' indicates that the specified region is in frame, i.e. that its first base corresponds to the first base of a codon. '1' indicates that there is one extra base, i.e. that the second base of the region corresponds to the first base of a codon, and '2' means that the third base of the region is the first base of a codon. If the strand is '-', then the first base of the region is value of <end>, because the corresponding coding region will run from <end> to <start> on the reverse strand.
    [group]
        An optional string-valued field that can be used as a name to group together a set of records. Typical uses might be to group the introns and exons in one gene prediction (or experimentally verified gene structure), or to group multiple regions of match to another sequence, such as an EST or a protein. See below for examples.
    
=> All strings (i.e. values of the <seqname>, <feature> or <group> fields) should be under 256 characters long, and should not include whitespace. The whole line should be under 32k long. A character limit is not very desirable, but helps write parsers in some languages. The slightly silly 32k limit is to allow plenty of space for comments/extra data.
=> Fields must be separated by TAB characters ('\t').

For a complete description of the format please refer to the link cited above.

Here an example of GFF1:  

    SEQ1	EMBL	atg	103	105	.	+	0
    SEQ1	EMBL	exon	103	172	.	+	0
    SEQ1	EMBL	splice5	172	173	.	+	.
    SEQ1	netgene	splice5	172	173	0.94	+	.
    SEQ1	genie	sp5-20	163	182	2.3	+	.
    SEQ1	genie	sp5-10	168	177	2.1	+	.
    SEQ2	grail	ATG	17	19	2.1	-   0
    SEQ2    pred

## GFF2 (become officialy the default version the 2000-9-29 but was proposed since 1998-12-16):
[Here a snapshot of the original page from SANGER (2000)](snapshots/sanger_gff2.md)
The GFF1 has evolved step by step and the 2000-9-29 the default version for GFF files became Version 2. Here we will see all the changes occured from the original version 1.

=> The **Gene Feature Finding** has been  generalized to accomodate to accommodate RNA and Protein feature files and has been renamed the **General Feature Format** while retaining the same acronym GFF.  

The main change from Version 1 to Version 2 is the addition of an optional 9th field with tag-value type structure (essentially semicolon-separated .ace format) used for any additional material on the line. Version 2 also allows '.' as a score, for features for which there is no score.
With the changes taking place to version 2 of the format, we also allow for feature sets to be defined over RNA and Protein sequences, as well as genomic DNA. This is used for example by the EMBOSS project to provide standard format output for all features as an option. In this case the <strand> and <frame> fields should be set to '.'. To assist this transition in specification, a new #Type Meta-Comment has been added.

## GTF (2002?)
Formats designed specifically for the human genome project. Created before June 2002 because it is mentioned in this paper: The Human Genome Browser at UCSC. Genome Res. 2002 Jun; 12(6): 996–1006. doi:  [10.1101/gr.229102]

## GTF2 (2003)
[Here the description from the Washington University in St. Louis](https://web.archive.org/web/20031212200757/http://genes.cse.wustl.edu/GTF2.html). Foudn from the Eval publication received the 18 july 2003 mentioning the address http://genes.cse.wustl.edu/GTF2.html that has been archived in the web-archive the 12/12/2003. Prior to the publication in BMV Bioinformatics (and after 1 January 2003 because it's the most recent  journal cited in his report) E. Kleiber released a Master project report named "Eval: A Gene Set Comparison System" where it mention and describe the GTF, maybe the first version of the format. 

## GTF2.2 (2007)
[Here the description from the Brent Lab (The Washington University in St. Louis) (http://mblab.wustl.edu/GTF22.html)

## GFF3

## Resume
