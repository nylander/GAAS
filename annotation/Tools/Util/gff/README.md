<h1 align="center">gff toolkit</h1>
---------------------------

Bench of tool to handle gff3 files.
To know more about gff3 format it's over there => https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

###################
# I) PREREQUISITE #
###################

Once you cloned the git repositpry you have to configure your installation as explained below:

## 1) Most of scripts use the BILS librairy which is located here: GAAS/annotation/BILS
Consequently, in order to use those scripts, you must add the library location to your path like that:

	export PERL5LIB=$PERL5LIB:/pathTo/GAAS/annotation

## 2) Bioperl must as well be installed as well as sepcific perl modules often not included by default (Moose,Clone)

### 2.1) Uppmax user (Swedish research cluster)=> This will be easy for you, It's a lucky day ;)
Just load the needed modules by executing these commands:<br>

	module load bioinfo-tools 
	module load BioPerl/1.6.922 
	module load perl_modules/5.18.4  # This module contains the Moose and Clone librairy

P.S: Instead to use the perl_modules/5.18.4 module you may also install your own local perl librairy by following these instructions: http://www.uppmax.uu.se/support/faq/software-faq/installing-local-perl-packages/  <br>

### 2.2) For non-uppmax user please follow these instructions:

During the below processes, if you encountered right problem, use the "super user" command to do the installation as administrator.

	e.g: sudo cpanm bioperl 
	
In that case you will be prompt to type the super user paswword. Obviously you need to have that privilege.

#### 2.2.1) Bioperl

The easiest way to have bioperl is to clone the bioperl-live project here => https://github.com/bioperl/bioperl-live
and add the path to it within your PERL5LIB path.

	export PERL5LIB=$PERL5LIB:/pathTo/bioperl-live 

Otherwise you should be able to install bioperl using your favorite package manager (cpan, cpanm, etc).

	e.g: cpanm bioperl

#### 2.2.2) Other mandatory modules
You must install other libraries (e.g Moose and Clone).
You can install them by using your  favorite package manager (cpan, cpanm, etc).

	e.g: cpanm Clone
	e.g: cpanm Moose 
	e.g: cpanm Graph::Directed
	e.g: cpanm LWP::UserAgent
	e.g: cpanm Statistics::R

#################################################
# II) Script name and classificaiton by prefix   #
					
_**As most as possible we will try to name the script with understandable names.
For that purpose we try to use a controled vocabulary**_

## A) gff vs gff3 prefix

Script not prefixed by gff3 but only with gff means that they havn't be checked or are not compatible with the gff3 standard. In other term, it means that a file not following the gff3 standards might not work with the script prefixed by gff3. Lot of modifcation could be post process if your file don't follow the gff3 standart. We will develop that in the part 3 of this readme.


## B) \_sq\_ AND \_sp\_

### B.1) \_sq\_ => Means SEQUENTIAL

The gff file is read and processed from top to the end. This is memory efficient !! 
But in other hand it hard to create complex script. Moreover, If data are not written sorted (e.g an exon of a gene located in the middle of the descritpion of another gene) some troubles could occur.

### B.2) \_sp\_ => Means SLURP

The gff file will be saved in memory before to process it. This is handle by the slurp_gff3_file_JD method. It has a memory cost. So if your gff3 files are over Gb size and your computer do not have enough ram memory, it might crash. 
That approach allows to peform more complicated task and more efficiency. Moreover, it allows to fix/correct, in the limit of the possibilities given by the format, the issues present in the gff you give in input. See part 3 for more information about it.


#################################################
# III) What does the SLURP method for you => GFF3 Standardization for a full GFF3 compliant to any tool !!!
#########
**_This method create a hash structure containing all the data in memory. We call it OMNISCIENT.<br>
The OMNISCNIENT structure is a three levels structure :_**

$omniscient{level1}{tag_l1}{level1_id} = feature <= tag could be gene,etc<br>
$omniscient{level2}{tag_l2}{idY} = @featureList <= tag could be mRNA,rRNA,tRNA,etc. idY is a level1_id (know as Parent attribute within the level2 feature). The @featureList is a list to be able to manage isoform cases.<br>
$omniscient{level3}{tag_l3}{idZ} =  @featureList <= tag could be exon,cds,utr3,utr5,etc. idZ is the ID of a level2 feature (know as Parent attribute within the level3 feature). The @featureList is a list to be able to put all the feature of a same tag together.<br>


It creates an ID attribute if missing <br>
It check for duplicated features (same position, same ID, same Parent)<br>
It expand level3 features (e.g. exon) sharing multiple mRNA (Parent attributes contains multiple parental mRNA). One exon by parental mRNA will be created.<br>
If a level 2 feature  doesn t have parent feature but has the attribute we create the level1 feature.<br>
If a feature  doesn t have the parent attribute we create the attribute !! But not the feature in the case of a parent of a level 3 feature. ( Could be implemented )<br>

INFO:Access to element of an omniscient is most of time from level1 to level3. Consequently if a level3 feature don't have any parent,  it will not be printed.<br>
