<h1 align="center">gff toolkit</h1>
---------------------------

Bench of tool to handle gff3 files.

###################
# 1) PREREQUISITE #
###################

1) Most of tool use the BILS librairy that can be found here: GAAS/annotation/BILS
Consequently, in order to use those script, you must add it to your path like that:

export PERL5LIB=$PERL5LIB:/pathTo/GAAS/annotation

2) Bioperl must as well be installed.

3) Some specific perl module like Clone.pm have to be installed too.

#################################################
# 2) Script name and classificaiton by prefix   #
#											    #
# As most as possible we will try to name the script with understandable names.
# For that purpose we try to use a controled vocabulary

A) Script not prefixed by gff3 but only with gff means that they havn't be checked or are not compatible with the gff3 standard. In other term, it means that a file not following the gff3 standards might not work with the script prefixed by gff3. Lot of modifcation could be post process if your file don't follow the gff3 standart. We will develop that in the part 3 of this readme.


B) _sq_ / _sp_
B.1) _sq_ => Means SEQUENTIAL = The gff file is read and processed from top to the end. This is memory efficient !! 
							 But in other hand it hard to create complex script. Moreover, If data are not written sorted (e.g an exon of a gene located in the middle of the descritpion of another gene) some troubles could occur.

B.2)_sp_ => Means SLURP = The gff file will be saved in memory before to process it. This is handle by the slurp_gff3_file_JD method. It has a memory cost. So if your gff3 files are over Gb size and your computer do not have enough ram memory, it might crash. 
That approach allows to peform more complicated task and more efficiency. Moreover, it allows to fix/correct, in the limit of the possibilities given by the format, the issues present in the gff you give in input. See part 3 for more information about it.


#################################################
# 3) What does the SLURP method for you
#########
This method create a hash structure containing all the data in memory. We call it OMNISCIENT
The OMNISCNIENT structure is a thre level structure :

$omniscient{level1}{level1_tag}{level1_id} = feature <= tag could be gene,etc
$omniscient{level2}{tagY}{idY} = @featureList <= tag could be mRNA,rRNA,tRNA,etc. idY is a level1_id (know as Parent attribute within the level2 feature). The @featureList is a list to be able to manage isoform cases.
$omniscient{level3}{tagZ}{idZ} =  @featureList <= tag could be exon,cds,utr3,utr5,etc. idZ is the ID of a level2 feature (know as Parent attribute within the level3 feature). The @featureList is a list to be able to put all the feature of a same tag together .


It creates an ID attribute if missing (It)
It check for duplicated features (same position, same ID, same Parent)
It expand level3 features (e.g. exon) sharing multiple mRNA (Parent attributes contains multiple parental mRNA). One exon by parental mRNA will be created.
If a level 2 feature  doesn t have parent feature but has the attribute we create the level1 feature.
If a feature  doesn t have the parent attribute we create the attribute !! But not the feature in the case of a parent of a level 3 feature. ( Ciould be implemented )

INFO:Access to element of an omniscient is most of time from level1 to level3. Consequently if a level3 feature don't have any parent,  it will not be printed.
