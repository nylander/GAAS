
# Mitochondrial annotation

List of tools that can be used for Mitochondrial annotation:


| year	| Tool name | Publication | Type	| Method | Organism | Comments | Output Format |
| --- | --- | --- | --- | --- | --- | --- | --- |
2011 | MAKER2 | Holt, C. & Yandell, M. MAKER2: an annotation pipeline and genome-database management tool for second-generation genome projects. BMC Bioinformatics 12, 491 (2011). | | | | 184 | It uses proteins, transcripts ... Abinitio: Augustus, Fgnesh,Genemark,snap. Need to modify the code to accept specific codon table for mitochondria. Not the best choice but provide nice protein alignments useful for manual curation |	
2013 | MitoAnnotator | Iwasaki W, Fukunaga T, Isagozawa R, Yamada K, Maeda, Y, Satoh TP, Sado T, Mabuchi K, Takeshima H, Miya M and Nishida M (2013) MitoFish and MitoAnnotator: a mitochondrial genome database of fish with an accurate and automatic annotation pipeline. Molecular Biology and Evolution 30: 2531-2540 | Ab initio (sensors + Neural network) | | Fish | MitoFish is a comprehensive and standardized fish mitochondrial genome database used by MitoAnnotator |
2014 | Prokka | Seemann T., Prokka: rapid prokaryotic genome annotation. Bioinformatics 2014 Jul 15;30(14):2068-9. PMID:24642063 | pipeline | Ab initio + evidence-based for functional annotation | prokaryote | | https://github.com/tseemann/prokka Do structural and functional annotation. No intron allowed! | .gff, .gbk, .fna, .faa, .ffn, .sqn, .fsa, .tbl, .err, .log, .txt, .tsv
| year	| Tool name | Publication | Type	| Method | Organism | Comments | Output Format |

Geseq which is web based tool and easy to use. It works well with both mitochondria and chloroplast.

AGORA - organelle genome annotation
Jung J, Kim JI, Jeong Y-S and Yi, G (2018) AGORA: organellar genome annotation from the amino acid and nucleotide references. Bioinformatics 34: 2661-2663

MFannot - organelle genome annotation
Developed by the labs of B. Franz Lang and Gertraud Burger (unpublished)
output not easy to deal with. Here the tool we use to convert in gff mfannot2gff.pl. But introns type II are not in that output.



MITOFY - plant mitochondrial genome annotation
Alverson AJ, Wei X, Rice DW, Stern DB, Barry K and Palmer JD (2010) Insights into the evolution of mitochondrial genome size from complete sequences of Citrullus lanatus and Cucurbita pepo (Cucurbitaceae). Molecular Biology and Evolution 27: 1436-1448

MITOS - metazoan mitochondrial genome annotation
Bernt M, Donath A, Jühling F, Externbrink F, Florentz C, Fritzsch G, Pütz J, Middendorf M and Stadler PF (2013) MITOS: improved de novo metazoan mitochondrial genome annotation. Molecular Phylogenetics and Evolution 69: 313-319

MOSAS - insect mitochondrial genome annotation
Sheffield NC, Hiatt KD, Valentine MC, Song H and Whiting MF (2010) Mitochondrial genomics in Orthoptera using MOSAS. Mitochondrial DNA 21: 87-104

(Nice resource where we catch most of those tools https://chlorobox.mpimp-golm.mpg.de/Alternative-Tools.html)
