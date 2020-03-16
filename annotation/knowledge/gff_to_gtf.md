# GFF to GTF conversion
## Review of the main conversion tools

It exists many GFF formats and many GTF formats 
(see [here](https://github.com/NBISweden/GAAS/blob/master/annotation/knowledge/gxf.md) for a complete review) and many tools
to perform the conversion. We will try to see in this review the main differences.

# Table of Contents

 * [The GFF file to convert](#the-gff-file-to-convert)
 * [The converters](#the-converters)
   * [AGAT](#agat-agat)
   * [gffread](#gffread) 
   * [GenomeTools](#genometools)
   * [ea-utils](#ea-utils)
   * [TransDecoder](#transdecoder)
   * [Kent utils](#kent-utils)
      
### The GFF file to convert

The test file is a GFF3 file:
```
##gff-version 3
##This is a test sample
scaffold625	maker	gene	337818	343277	.	+	.	ID=CLUHARG00000005458;Name=TUBB3_2
scaffold625	maker	mRNA	337818	343277	.	+	.	ID=CLUHART00000008717;Parent=CLUHARG00000005458
scaffold625	maker	tss	337915	337918	.	+	.	ID=CLUHART00000008717:tss;Parent=CLUHART00000008717
scaffold625	maker	CDS	337915	337971	.	+	0	ID=CLUHART00000008717:cds;Parent=CLUHART00000008717
scaffold625	maker	CDS	340733	340841	.	+	0	ID=CLUHART00000008717:cds;Parent=CLUHART00000008717
scaffold625	maker	CDS	341518	341628	.	+	2	ID=CLUHART00000008717:cds;Parent=CLUHART00000008717
scaffold625	maker	CDS	341964	343033	.	+	2	ID=CLUHART00000008717:cds;Parent=CLUHART00000008717
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

### ea-utils

### TransDecoder

### Kent utils
