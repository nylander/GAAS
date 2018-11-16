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

## GFF1(1997):
[Here a snapshot of the olderst description of the format I found (2000)]((http://htmlpreview.github.io/?https://github.com/NBISweden/GAAS/blob/master/annotation/CheatSheet/snapshots/GFF_Spec.html).

This GFF1 format was containing 8 fields separated by tabulations. Fields are:  

    <seqname> <source> <feature> <start> <end> <score> <strand> <frame>

Here an example of GFF1:  

    SEQ1	EMBL	atg	103	105	.	+	0
    SEQ1	EMBL	exon	103	172	.	+	0
    SEQ1	EMBL	splice5	172	173	.	+	.
    SEQ1	netgene	splice5	172	173	0.94	+	.
    SEQ1	genie	sp5-20	163	182	2.3	+	.
    SEQ1	genie	sp5-10	168	177	2.1	+	.
    SEQ2	grail	ATG	17	19	2.1	- 0

## GFF2 (become officialy the default version the 2000-9-29 but was proposed since 1998-12-16):
[Here a snapshot of the original page from SANGER (2000)](snapshots/sanger_gff2.md)
2000-9-29 The default version for GFF files is now Version 2. The **Gene Feature Finding** has been  generalized to accomodate to accommodate RNA and Protein feature files and has been renamed the **General Feature Format** while retaining the same acronym GFF.  

The main change from Version 1 to Version 2 is the requirement for a tag-value type structure (essentially semicolon-separated .ace format) for any additional material on the line, following the mandatory fields. Version 2 also allows '.' as a score, for features for which there is no score.

## GTF (2002?)
Formats designed specifically for the human genome project. Created before June 2002 because it is mentioned in this paper: The Human Genome Browser at UCSC. Genome Res. 2002 Jun; 12(6): 996–1006. doi:  [10.1101/gr.229102]

## GTF2 (2003)
[Here the description from the Washington University in St. Louis](https://web.archive.org/web/20031212200757/http://genes.cse.wustl.edu/GTF2.html). Foudn from the Eval publication received the 18 july 2003 mentioning the address http://genes.cse.wustl.edu/GTF2.html that has been archived in the web-archive the 12/12/2003. Prior to the publication in BMV Bioinformatics (and after 1 January 2003 because it's the most recent  journal cited in his report) E. Kleiber released a Master project report named "Eval: A Gene Set Comparison System" where it mention and describe the GTF, maybe the first version of the format. 

## GTF2.2 (2007)
[Here the description from the Brent Lab (The Washington University in St. Louis) (http://mblab.wustl.edu/GTF22.html)
