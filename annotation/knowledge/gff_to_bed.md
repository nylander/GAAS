# GFF to BED conversion
## Review of the main conversion tools

It exists many GFF formats and many GTF formats 
(see [here](https://github.com/NBISweden/GAAS/blob/master/annotation/knowledge/gxf.md) for a complete review) and many tools
to perform the conversion. We will try to see in this review the main differences.

# Table of Contents

 * [Test resume](#test-resume)
 * [The GFF file to convert](#the-gff-file-to-convert)
 * [The converters](#the-converters)
   * [AGAT](#agat)
   * [gffread](#gffread) 
   * [GenomeTools](#genometools)
   * [ea-utils](#ea-utils)
   * [TransDecoder](#transdecoder)
   * [Kent utils](#kent-utils)
 * [Feature types in GTF versions](feature-types-in-gtf-versions)

### Test resume

tool | respect GTF format | UTR conserved | attribute conserved | Comment
-- | -- | -- | -- | -- |
[AGAT](https://github.com/NBISweden/AGAT) | Yes - All (default GTF3) | Yes it converts UTR terms to the appropriate ones according to the GTF version selected.| Yes - All | Can take any GTF GFF as input. The only one keeping comments at the beginning of the file.
[gffread](https://github.com/gpertea/gffread) | No - They say GTF2.2 but it is not: transcript should be removed; start_codon and stop_codon should stay. | No | No  | 
[GenomeTools](https://github.com/genometools/genometools) | No - only CDS and exon kept | No | No | gene_id and transcript_id get new identifiers.
[ea-utils](https://github.com/ExpressionAnalysis/ea-utils) |  No - only CDS and exon kept | No | No |
[TransDecoder](https://github.com/TransDecoder/TransDecoder) |  No - start and stop codon removed | No | Yes - Name only | Needs the fasta file for the conversion. 
[Kent utils](http://hgdownload.cse.ucsc.edu/admin/exe/) | No - gene is missing or transcript is superfluous to be compliant to one of the GTF format | No | No | Create a new attribute 'gene_name'

### The GFF file to convert

The test file is a GFF3 file:

```
##gff-version 3
##This is a test sample
scaffold625	maker	gene	337818	343277	.	+	.	ID=CLUHARG00000005458;Name=TUBB3_2
scaffold625	maker	mRNA	337818	343277	.	+	.	ID=CLUHART00000008717;Parent=CLUHARG00000005458
scaffold625	maker	tss	337916	337918	.	+	.	ID=CLUHART00000008717:tss;Parent=CLUHART00000008717
scaffold625	maker	start_codon	337916	337918	.	+	.	ID=CLUHART00000008717:start;Parent=CLUHART00000008717
scaffold625	maker	CDS	337915	337971	.	+	0	ID=CLUHART00000008717:cds;Parent=CLUHART00000008717
scaffold625	maker	CDS	340733	340841	.	+	0	ID=CLUHART00000008717:cds;Parent=CLUHART00000008717
scaffold625	maker	CDS	341518	341628	.	+	2	ID=CLUHART00000008717:cds;Parent=CLUHART00000008717
scaffold625	maker	CDS	341964	343033	.	+	2	ID=CLUHART00000008717:cds;Parent=CLUHART00000008717
scaffold625	maker	stop_codon	343031	343033	.	+	.	ID=CLUHART00000008717:stop;Parent=CLUHART00000008717
scaffold625	maker	exon	337818	337971	.	+	.	ID=CLUHART00000008717:exon1;Parent=CLUHART00000008717
scaffold625	maker	exon	340733	340841	.	+	.	ID=CLUHART00000008717:exon2;Parent=CLUHART00000008717
scaffold625	maker	exon	341518	341628	.	+	.	ID=CLUHART00000008717:exon3;Parent=CLUHART00000008717
scaffold625	maker	exon	341964	343277	.	+	.	ID=CLUHART00000008717:exon4;Parent=CLUHART00000008717
scaffold625	maker	five_prime_utr	337818	337914	.	+	.	ID=CLUHART00000008717:five_prime_utr;Parent=CLUHART00000008717
scaffold625	maker	three_prime_UTR	343034	343277	.	+	.	ID=CLUHART00000008717:three_prime_utr;Parent=CLUHART00000008717
```

### AGAT

AGAT v0.2.2  

`agat_convert_sp_gff2gtf.pl --gff 1_test.gff -o 1_test_agat.gtf`

```
##gtf-version 3
##This is a test sample
scaffold625	maker	gene	337818	343277	.	+	.	ID CLUHARG00000005458; Name TUBB3_2; gene_id CLUHARG00000005458
scaffold625	maker	transcript	337818	343277	.	+	.	ID CLUHART00000008717; Parent CLUHARG00000005458; gene_id CLUHARG00000005458; original_biotype mrna; transcript_id CLUHART00000008717
scaffold625	maker	exon	337818	337971	.	+	.	ID "CLUHART00000008717:exon1"; Parent CLUHART00000008717; gene_id CLUHARG00000005458; transcript_id CLUHART00000008717
scaffold625	maker	exon	340733	340841	.	+	.	ID "CLUHART00000008717:exon2"; Parent CLUHART00000008717; gene_id CLUHARG00000005458; transcript_id CLUHART00000008717
scaffold625	maker	exon	341518	341628	.	+	.	ID "CLUHART00000008717:exon3"; Parent CLUHART00000008717; gene_id CLUHARG00000005458; transcript_id CLUHART00000008717
scaffold625	maker	exon	341964	343277	.	+	.	ID "CLUHART00000008717:exon4"; Parent CLUHART00000008717; gene_id CLUHARG00000005458; transcript_id CLUHART00000008717
scaffold625	maker	CDS	337915	337971	.	+	0	ID "CLUHART00000008717:cds"; Parent CLUHART00000008717; gene_id CLUHARG00000005458; transcript_id CLUHART00000008717
scaffold625	maker	CDS	340733	340841	.	+	0	ID "CLUHART00000008717:cds"; Parent CLUHART00000008717; gene_id CLUHARG00000005458; transcript_id CLUHART00000008717
scaffold625	maker	CDS	341518	341628	.	+	2	ID "CLUHART00000008717:cds"; Parent CLUHART00000008717; gene_id CLUHARG00000005458; transcript_id CLUHART00000008717
scaffold625	maker	CDS	341964	343033	.	+	2	ID "CLUHART00000008717:cds"; Parent CLUHART00000008717; gene_id CLUHARG00000005458; transcript_id CLUHART00000008717
scaffold625	maker	five_prime_utr	337818	337914	.	+	.	ID "CLUHART00000008717:five_prime_utr"; Parent CLUHART00000008717; gene_id CLUHARG00000005458; transcript_id CLUHART00000008717
scaffold625	maker	start_codon	337916	337918	.	+	.	ID "CLUHART00000008717:start"; Parent CLUHART00000008717; gene_id CLUHARG00000005458; transcript_id CLUHART00000008717
scaffold625	maker	stop_codon	343031	343033	.	+	.	ID "CLUHART00000008717:stop"; Parent CLUHART00000008717; gene_id CLUHARG00000005458; transcript_id CLUHART00000008717
scaffold625	maker	three_prime_utr	343034	343277	.	+	.	ID "CLUHART00000008717:three_prime_utr"; Parent CLUHART00000008717; gene_id CLUHARG00000005458; original_biotype three_prime_UTR; transcript_id CLUHART00000008717
```

### gffread

gffread 0.11.4

`gffread -E 1_test.gff -T -o  1_test_gffread.gtf` 

```
scaffold625	maker	transcript	337818	343277	.	+	.	transcript_id "CLUHART00000008717"; gene_id "CLUHARG00000005458";
scaffold625	maker	exon	337818	337971	.	+	.	transcript_id "CLUHART00000008717"; gene_id "CLUHARG00000005458";
scaffold625	maker	exon	340733	340841	.	+	.	transcript_id "CLUHART00000008717"; gene_id "CLUHARG00000005458";
scaffold625	maker	exon	341518	341628	.	+	.	transcript_id "CLUHART00000008717"; gene_id "CLUHARG00000005458";
scaffold625	maker	exon	341964	343277	.	+	.	transcript_id "CLUHART00000008717"; gene_id "CLUHARG00000005458";
scaffold625	maker	CDS	337915	337971	.	+	0	transcript_id "CLUHART00000008717"; gene_id "CLUHARG00000005458";
scaffold625	maker	CDS	340733	340841	.	+	0	transcript_id "CLUHART00000008717"; gene_id "CLUHARG00000005458";
scaffold625	maker	CDS	341518	341628	.	+	2	transcript_id "CLUHART00000008717"; gene_id "CLUHARG00000005458";
scaffold625	maker	CDS	341964	343033	.	+	2	transcript_id "CLUHART00000008717"; gene_id "CLUHARG00000005458";
```

### GenomeTools

GenomeTools 1.6.1  
The help says it convert into GTF2.2

`gt gff3_to_gtf 1_test.gff > 1_test_genometools.gtf`

```
scaffold625	maker	exon	337818	337971	.	+	.	gene_id "1"; transcript_id "1.1";
scaffold625	maker	exon	340733	340841	.	+	.	gene_id "1"; transcript_id "1.1";
scaffold625	maker	exon	341518	341628	.	+	.	gene_id "1"; transcript_id "1.1";
scaffold625	maker	exon	341964	343277	.	+	.	gene_id "1"; transcript_id "1.1";
scaffold625	maker	CDS	337915	337971	.	+	0	gene_id "1"; transcript_id "1.1";
scaffold625	maker	CDS	340733	340841	.	+	0	gene_id "1"; transcript_id "1.1";
scaffold625	maker	CDS	341518	341628	.	+	2	gene_id "1"; transcript_id "1.1";
scaffold625	maker	CDS	341964	343033	.	+	2	gene_id "1"; transcript_id "1.1";
```

### ea-utils

[ea-utils](https://github.com/ExpressionAnalysis/ea-utils) commit 2b3d8c5d148801c98a2b3f3d54009a72c5b99521

`./gff2gtf-eautils test_1.gff >  1_test_ea-utils.gtf`

```
scaffold625	maker	exon	337818	337971	0	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717:CLUHARG00000005458";
scaffold625	maker	CDS	337915	337971	0	+	0	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717:CLUHARG00000005458";
scaffold625	maker	CDS	340733	340841	0	+	0	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717:CLUHARG00000005458";
scaffold625	maker	exon	340733	340841	0	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717:CLUHARG00000005458";
scaffold625	maker	CDS	341518	341628	0	+	2	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717:CLUHARG00000005458";
scaffold625	maker	exon	341518	341628	0	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717:CLUHARG00000005458";
scaffold625	maker	CDS	341964	343033	0	+	2	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717:CLUHARG00000005458";
scaffold625	maker	exon	341964	343277	0	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717:CLUHARG00000005458";
```

### TransDecoder

Transdecoder  v5.5.0

`gff3_gene_to_gtf_format.pl test_1.gff test_1.fa > 1_test_transdecoder.gtf`

```
scaffold625	maker	gene	337818	343277	0	+	.	gene_id "CLUHARG00000005458"; Name "TUBB3_2";
scaffold625	maker	transcript	337818	343277	0	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; Name "TUBB3_2";
scaffold625	maker	exon	337818	337971	0	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; Name "TUBB3_2";
scaffold625	maker	CDS	337818	337971	0	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; Name "TUBB3_2";
scaffold625	maker	exon	340733	340841	0	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; Name "TUBB3_2";
scaffold625	maker	CDS	340733	340841	0	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; Name "TUBB3_2";
scaffold625	maker	exon	341518	341628	0	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; Name "TUBB3_2";
scaffold625	maker	CDS	341518	341628	0	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; Name "TUBB3_2";
scaffold625	maker	exon	341964	343277	0	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; Name "TUBB3_2";
scaffold625	maker	CDS	341964	343277	0	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; Name "TUBB3_2";
```

### Kent utils

version from 26-Feb-2020

`./gff3ToGenePred.dms 1_test.gff temp.genePred`
`./genePredToGtf.dms file temp.genePred 1_test_genePred.gtf`

```
scaffold625	temp.genePred	transcript	337818	343277	.	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717";  gene_name "CLUHARG00000005458";
scaffold625	temp.genePred	exon	337818	337971	.	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; exon_number "1"; exon_id "CLUHART00000008717.1"; gene_name "CLUHARG00000005458";
scaffold625	temp.genePred	CDS	337915	337971	.	+	0	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; exon_number "1"; exon_id "CLUHART00000008717.1"; gene_name "CLUHARG00000005458";
scaffold625	temp.genePred	exon	340733	340841	.	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; exon_number "2"; exon_id "CLUHART00000008717.2"; gene_name "CLUHARG00000005458";
scaffold625	temp.genePred	CDS	340733	340841	.	+	0	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; exon_number "2"; exon_id "CLUHART00000008717.2"; gene_name "CLUHARG00000005458";
scaffold625	temp.genePred	exon	341518	341628	.	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; exon_number "3"; exon_id "CLUHART00000008717.3"; gene_name "CLUHARG00000005458";
scaffold625	temp.genePred	CDS	341518	341628	.	+	2	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; exon_number "3"; exon_id "CLUHART00000008717.3"; gene_name "CLUHARG00000005458";
scaffold625	temp.genePred	exon	341964	343277	.	+	.	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; exon_number "4"; exon_id "CLUHART00000008717.4"; gene_name "CLUHARG00000005458";
scaffold625	temp.genePred	CDS	341964	343030	.	+	2	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; exon_number "4"; exon_id "CLUHART00000008717.4"; gene_name "CLUHARG00000005458";
scaffold625	temp.genePred	start_codon	337915	337917	.	+	0	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; exon_number "1"; exon_id "CLUHART00000008717.1"; gene_name "CLUHARG00000005458";
scaffold625	temp.genePred	stop_codon	343031	343033	.	+	0	gene_id "CLUHARG00000005458"; transcript_id "CLUHART00000008717"; exon_number "4"; exon_id "CLUHART00000008717.4"; gene_name "CLUHARG00000005458";
```

# The bed format

Detailed information can be found here: [https://genome.ucsc.edu/FAQ/FAQformat.html](https://genome.ucsc.edu/FAQ/FAQformat.html)  
Below a description of the different fields:

column | feature type | mandatory | comment
-- | -- | -- | -- |
1 | chrom | X |  The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
2 | chromStart | X | The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
3 | chromEnd | X | The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
4 | name  |  | Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
5 | score |  | A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray).
6 | strand |  | Defines the strand - either '+' or '-'.
7 | thickStart |  | The starting position at which the feature is drawn thickly
8 | thickEnd |  | The ending position at which the feature is drawn thickly
9 | itemRgb |  | An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.
10 | blockCount |  | The number of blocks (exons) in the BED line.
11 | blockSizes |  | A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
12 | blockStarts |  | A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount.


/!\ location BED format is 0-based, half-open [start-1, end), while GFF is 1-based, closed [start, end].
