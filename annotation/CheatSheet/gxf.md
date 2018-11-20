The gene-finding format / general feature format  
GFF / GFF1 / GFF2 / GFF2.5 / GFF3 / GTF / GTF2.2 formats
===========================

It's often hard to understand and differentiate all GFF/GTF formats/flavors. Here is an overview of the format and its history to help disentangle this complexity.


## Forewords:
⇨	When I use the term gff it includes all gff formats/flavors. (The first version of the format was not called gff1 but gff. But to make it easier I will always express the version of the format saying gff1 when it's the first version of it. So from now when I say gff it means all gff formats/flavors).  
⇨	I the same way, when I use the term gtf it includes all gtf formats/flavors.  
⇨	I have created the term **gxf** that means all the gff and gtf formats/flavors.

The GFF Protocol Specification was initially proposed by **Richard Durbin** and **David Haussler**.
The format has been originaly developed to help the gene prediction (or gene finding) world. Indeed the gene prediction methods are based on two main steps, **first** finding signals (starts, splice sites, stops, motifs, etc.) and regions (exons, introns, protein domains etc.); **secondly** combining these to give complete gene, RNA transcript or protein structures. These two steps were usually performed within the same program. In order to decoupled them they have created the format called GFF ('Gene-Finding Format') allowing the transfer of feature information from a tool to another one.

The GFF fomat has been developed to be easy to parse and process by a variety of programs in different languages (e.g Unix tools as grep and sort, perl, awk, etc). For these reasons, they decided that each feature is described on a single line, and line order is not relevant.

A GFF record is an extension of a basic (name,start,end) tuple (or "NSE") that can be used to identify a substring of a biological sequence. 

## GFF0(before 13-11-1997)

There is no clear information about how look the format at that time but it was close to the GFF1 format specification without the field "source" added the 1997-11-13.

## GFF1(13-11-1997):

For a complete description of the format please refer to this link:
[https://web.archive.org/web/19980222142332/http://www.sanger.ac.uk:80/~rd/gff.html](https://web.archive.org/web/19980222142332/http://www.sanger.ac.uk:80/~rd/gff.html). This is the oldest description of the format I found (1998-02-22).

I consider the format as GFF1 when they definitly defined the 9th field of the format (1997-11-13 rd: added extra "source" field as discussed at Newton Institute meeting 971029). Before that the format was existing but was at the stage of version 0.

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

Extra features of the format:

    Comments

    Comments are allowed, starting with "#" as in Perl, awk etc. Everything following # until the end of the line is ignored. Effectively this can be used in two ways. Either it must be at the beginning of the line (after any whitespace), to make the whole line a comment, or the comment could come after all the required fields on the line.
    We also permit extra information to be given on the line following the group field without a '#' character. This allows extra method-specific information to be transferred with the line. However, we discourage overuse of this feature: better to find a way to do it with more true feature lines, and perhaps groups.

    ## comment lines for meta information

    There is a set of standardised (i.e. parsable) ## line types that can be used optionally at the top of a gff file. The philosophy is a little like the special set of %% lines at the top of postscript files, used for example to give the BoundingBox for EPS files.
    Current proposed ## lines are:

     ##gff-version 1 
    GFF version - in case it is a real success and we want to change it. The current version is 1.
     ##source-version {source} {version text} 
    So that people can record what version of a program or package was used to make the data in this file. I suggest the version is text without whitespace. That allows things like 1.3, 4a etc.
     ##date {date} 
    The date the file was made, or perhaps that the prediction programs were run. We suggest to use astronomical format: 1997-11-08 for 8th November 1997, first because these sort properly, and second to avoid any US/European bias.

     ##DNA {seqname}
     ##acggctcggattggcgctggatgatagatcagacgac
     ##...
     ##end-DNA
    To give a DNA sequence. Several people have pointed out that it may be convenient to include the sequence in the file. It should not become mandatory to do so. Often the seqname will be a well-known identifier, and the sequence can easily be retrieved from a database, or an accompanying file.
     ##sequence-region {seqname} {start} {end} 
    To indicate that this file only contains entries for the the specified subregion of a sequence.
    Please feel free to propose new ## lines. The ## line proposal came out of some discussions including Anders Krogh, David Haussler, people at the Newton Institute on 1997-10-29 and some email from Suzanna Lewis. Of course, naive programs can ignore all of these...

Here an example of GFF1:  

    ##gff-version 1 
    SEQ1	EMBL	atg	103	105	.	+	0
    SEQ1	EMBL	exon	103	172	.	+	0
    SEQ1	EMBL	splice5	172	173	.	+	.
    SEQ1	netgene	splice5	172	173	0.94	+	.
    # this is comment that will be skipped by the parser
    SEQ1	genie	sp5-20	163	182	2.3	+	.
    SEQ1	genie	sp5-10	168	177	2.1	+	.
    SEQ2	grail	ATG	17	19	2.1	-   0
    SEQ3    pred    exon    100 135 .   +   0   locus1 # this is also a comment that will be skipped by the parser
    SEQ3    pred    exon    235 260 .   +   2   locus1 This is an example of extra information... They discourage overuse of this feature.
    SEQ3    pred    exon    360 396 .   +   0   locus1

## GFF2 (29-09-2000):
/!\ Note: Some of the changes we will see have been implemented before the offical release of the version 2. As consequence, several interemediate states between version 1 and 2 have existed. We can call them GFF1.X. I will not discuss further these intermediate states.

**16/12/98**: Discussions with **Lincoln Stein** and **others**,the Version 2 format of GFF is proposed.  
**17/11/99**: 'Gene Feature Finding' Version 2 format is conceptually generalized to be the 'General Feature Format'

The GFF2 format is conceptualized since the 16/12/98 but becomes officially the default version the 2000-9-29.
[Here is the official description](snapshots/sanger_gff2.md) which is a snapshot from here: https://web.archive.org/web/20010208224442/http://www.sanger.ac.uk:80/Software/formats/GFF/GFF_Spec.shtml.  
You can find the first description (03 Feb 2000) of the GFF2 [here] (snapshots/GFF2_Spec_first_draft_03_feb_2000.html) that comes from here:   ftp://ftp.sanger.ac.uk/pub/resources/software/gff-old/gff/).

Here we will review changes occured from the version 1. 

=> The **Gene Feature Finding** has been  generalized to accomodate to accommodate RNA and Protein feature files and has been renamed the **General Feature Format** while retaining the same acronym GFF.  

The main change from Version 1 to Version 2 is the addition of an optional 9th field with tag-value type structure (essentially semicolon-separated .ace format) used for any additional material on the line. Version 2 also allows '.' as a score, for features for which there is no score.
With the changes taking place to version 2 of the format, we also allow for feature sets to be defined over RNA and Protein sequences, as well as genomic DNA. This is used for example by the EMBOSS project to provide standard format output for all features as an option. In this case the <strand> and <frame> fields should be set to '.'. To assist this transition in specification, a new #Type Meta-Comment has been added.

Definition

This GFF2 format contains 8 madatory fields and 9th one optional. Fields are:  

      <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [group/attributes] [comments]

Definition of these fields are (For better lisibility here is included only differences with GFF1):

    <seqname>
        /
    <source>
        /
    <feature>
        Version 2 change: Standard Table of Features - we would like to enforce a standard nomenclature for common GFF features. This does not forbid the use of other features, rather, just that if the feature is obviously described in the standard list, that the standard label should be used. For this standard table we propose to fall back on the international public standards for genomic database feature annotation, specifically, the DDBJ/EMBL/GenBank feature table.
    <start>, <end>
        Version 2 change: version 2 condones values of <start> and <end> that extend outside the reference sequence. This is often more natural when dumping from acedb, rather than clipping. It means that some software using the files may need to clip for itself.
    <score>
        Version 2 change: When there is no score (i.e. for a sensor that just records the possible presence of a signal, as for the EMBL features above) you should use '.' instead of 0.
    <strand>
        Version 2 change: This field is left empty '.' for RNA and protein features.
    <frame>
        Version 2 change: This field is left empty '.' for RNA and protein features.
    [group/attribute]
        [New] Standard Table of Attribute Tag Identifiers The semantics of tags in attribute field tag-values pairs has not yet been completely formalized, however a useful constraint is that they be equivalent, where appropriate, to DDBJ/EMBL/GenBank feature 'qualifiers' of given features (see EMBL feature descriptions).

    In addition to these, ACEDB typically dumps GFF with specific tag-value pairs for given feature types. These tag-value pairs may be considered 'standard' GFF tag-values with respect to ACEDB databases. (rbsk: These will be summarized in a table here in the near future)

    Version 2 change: In version 2, the optional [group] field is renamed to [attribute] (09/99) and must have an tag value structure following the syntax used within objects in a .ace file, flattened onto one line by semicolon separators. Tags must be standard identifiers ([A-Za-z][A-Za-z0-9_]*). Free text values must be quoted with double quotes. Note: all non-printing characters in such free text value strings (e.g. newlines, tabs, control characters, etc) must be explicitly represented by their C (UNIX) style backslash-escaped representation (e.g. newlines as '\n', tabs as '\t'). As in ACEDB, multiple values can follow a specific tag. The aim is to establish consistent use of particular tags, corresponding to an underlying implied ACEDB model if you want to think that way (but acedb is not required). Examples of these would be:
    seq1     BLASTX  similarity   101  235 87.1 + 0	Target "HBA_HUMAN" 11 55 ; E_value 0.0003
    dJ102G20 GD_mRNA coding_exon 7105 7201   .  - 2 Sequence "dJ102G20.C1.1"

    => Version 2 change: field and line size limitations are removed; however, fields (except the optional [attribute] field above) must still not include whitespace.
    => Version 2 note: previous Version 2 permission to use arbitrary whitespace as field delimiters is now revoked! (99/02/26)

Extra features of the format:

    Comments

    [...]
    We also permit extra information to be given on the line following the attribute field without a '#' character (Version 2 change: this extra information must be delimited by the '#' comment delimiter OR by another tab field delimiter character, following any and all [attribute] field tag-value pairs).
    [...]

    Version 2 change: we gave in and defined a structured way of passing additional information, as described above under [attribute]. But the sentiment of this paragraph still applies - don't overuse the tag-value syntax. The use of tag-value pairs (with whitespace) renders problematic the parsing of Version 1 style comments (following the attribute field, without a '#' character), so in Version 2, such [attribute] trailing comments must either start with the "#" as noted above, or with at least one additional tab character. Moreover, '#' characters embedded within quoted text string values of [attribute] tag-values should not be parsed as the beginning of a comment.

    ## comment lines for meta information

    There is a set of standardised (i.e. parsable) ## line types that can be used optionally at the top of a gff file. The philosophy is a little like the special set of %% lines at the top of postscript files, used for example to give the BoundingBox for EPS files.
    Current proposed ## lines are:

      ##gff-version 2 
        The current version is 2. (Version 2 change!)
     ##source-version {source} {version text} 
        /
     ##date {date} 
        /
      ##Type <type> [<name>] 
        [New] The type of host sequence described by the features. Standard types are 'DNA', 'Protein' and 'RNA'. The optional <name> allows multiple ##Type definitions describing multiple GFF sets in one file, each which have a distinct type. If the name is not provided, then all the features in the file are of the given type. Thus, with this meta-comment, a single file could contain DNA, RNA and Protein features, for example, representing a single genomic locus or 'gene', alongside type-specific features of its transcribed mRNA and translated protein sequences. If no ##Type meta-comment is provided for a given GFF file, then the type is assumed to be DNA.

     ##DNA {seqname}
     ##acggctcggattggcgctggatgatagatcagacgac
     ##...
     ##end-DNA
        /
        
     ##RNA <seqname>
     ##acggcucggauuggcgcuggaugauagaucagacgac
     ##...
     ##end-RNA
        Similar to DNA. Creates an implicit ##Type RNA <seqname> directive.

     ##Protein <seqname>
     ##MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSF
     ##...
     ##end-Protein
        Similar to DNA. Creates an implicit ##Type Protein <seqname> directive.

      ##sequence-region {seqname} {start} {end} 
        /

Here an example of GFF2:  

    ##gff-version 2 
    SEQ1	EMBL	atg	103	105	.	+	0
    SEQ1	netgene	splice5	172	173	0.94	+	.
    # this is comment that will be skipped by the parser
    SEQ1	genie	sp5-10	168	177	2.1	+	.
    SEQ2	grail	ATG	17	19	2.1	-   0
    SEQ3    BLASTX    similarity    100 135 .   +   0   Target "HBA_HUMAN" ; E_value 0.0003 # this is also a comment that will be skipped by the parser
    SEQ3    BLASTX    similarity    235 260 .   +   2   Target "HBA_HUMAN" ; E_value 0.0005
    SEQ3    BLASTX    similarity    360 396 .   +   0   Target "HBA_HUMAN" ; E_value 0.001

## GTF (2002?)
Formats designed specifically for the human genome project. Created before June 2002 because it is mentioned in this paper: The Human Genome Browser at UCSC. Genome Res. 2002 Jun; 12(6): 996–1006. doi:  [10.1101/gr.229102]

## GTF2 (2003)
[Here the description from the Washington University in St. Louis](https://web.archive.org/web/20031212200757/http://genes.cse.wustl.edu/GTF2.html). Foudn from the Eval publication received the 18 july 2003 mentioning the address http://genes.cse.wustl.edu/GTF2.html that has been archived in the web-archive the 12/12/2003. Prior to the publication in BMV Bioinformatics (and after 1 January 2003 because it's the most recent  journal cited in his report) E. Kleiber released a Master project report named "Eval: A Gene Set Comparison System" where it mention and describe the GTF, maybe the first version of the format. 

## GTF2.2 (2007)
[Here the description from the Brent Lab (The Washington University in St. Louis) (http://mblab.wustl.edu/GTF22.html)

## GFF3 ()

GFF3 addresses several shortcomings in its predecessor GFF2.

## Resume
