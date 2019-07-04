
# Chloroplast annotation

List of tools that can be used

| year	| Tool name | Publication | Type	| Method | Organism | Comments | Output Format |
| --- | --- | --- | --- | --- | --- | --- | --- |
2004 | DOGMA | Wyman SK, Jansen RK and Boore JL (2004) Automatic annotation of organellar genomes with DOGMA. Bioinformatics 20: 3252-3255 | | | plant chloroplast and animal mitochondrial | |
2011 | MAKER2 | Holt, C. & Yandell, M. MAKER2: an annotation pipeline and genome-database management tool for second-generation genome projects. BMC Bioinformatics 12, 491 (2011). | | | | It uses proteins, transcripts ... Abinitio: Augustus, Fgnesh,Genemark,snap. Need to modify the code to accept specific codon table for mitochondria. Not the best choice but provide nice protein alignments useful for manual curation |	|
2012 | CpGAVAS | Liu C, Shi L, Zhu Y, Chen H, Zhang J, Lin X and Guan X (2012) CpGAVAS, an integrated web server for the annotation, visualization, analysis, and GenBank submission of completely sequenced chloroplast genome sequences. BMC Genomics 13: 715 | | | | web based tool |
2013 | CGAP | Cheng J, Zeng Xu, Ren G and Liu Z (2013) CGAP: a new comprehensive platform for the comparative analysis of chloroplast genomes. BMC Bioinformatics 14:95 | | | | |
2014 | Prokka | Seemann T., Prokka: rapid prokaryotic genome annotation. Bioinformatics 2014 Jul 15;30(14):2068-9. PMID:24642063 | pipeline | Ab initio + evidence-based for functional annotation | prokaryote | https://github.com/tseemann/prokka Do structural and functional annotation. No intron allowed! | .gff, .gbk, .fna, .faa, .ffn, .sqn, .fsa, .tbl, .err, .log, .txt, .tsv
2015 | ORG.Annot | Developed by Eric Coissac (unpublished) | Pipeline | | | |
2017 | Geseq | Tillich M, Lehwark P, Pellizzer T, Ulbricht-Jones ES, Fischer A, Bock R and Greiner S (2017) GeSeq – versatile and accurate annotation of organelle genomes. Nucleic Acids Research 45: W6-W11 | | | chloroplast mitochondria | web based tool | |
2019 | CpGAVAS2 | Shi L, Chen H, Jiang M, Wang L, Wu X, Huang L and Liu C (2019) CPGAVAS2, an integrated plastome sequence annotator and analyzer. Nucleic Acids Research, gkz345 | | | |
| year	| Tool name | Publication | Type	| Method | Organism | Comments | Output Format |



AGORA - organelle genome annotation
Jung J, Kim JI, Jeong Y-S and Yi, G (2018) AGORA: organellar genome annotation from the amino acid and nucleotide references. Bioinformatics 34: 2661-2663

MFannot - organelle genome annotation
Developed by the labs of B. Franz Lang and Gertraud Burger (unpublished)
output not easy to deal with. Here the tool we use to convert in gff mfannot2gff.pl. But introns type II are not in that output.


Plann - chloroplast genome annotation
Huang DI and Cronk QCB (2015) Plann: a command-line application for annotating plastome sequences (2015). Applications in Plant Sciences 3: 1500026

Verdant - chloroplast genome annotation and phylogenetic analysis
McKain MR, Hartsock RH, Wohl MM and Kellogg EA (2017) Verdant: automated annotation, alignment and phylogenetic analysis of whole chloroplast genomes. Bioinformatics 33: 130-132


(Nice resource where we catch most of those tools https://chlorobox.mpimp-golm.mpg.de/Alternative-Tools.html)
