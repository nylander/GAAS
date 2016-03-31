#!/usr/bin/env perl

###########################################
# gff manager v8 - Jacques Dainat 11/2014 # 
###########################################

#libraries
use strict;
use warnings;
use Data::Dumper;
use Carp;
use POSIX qw(strftime);
use Getopt::Long;
use IO::File;
use Pod::Usage;
use List::MoreUtils qw(uniq);
use Bio::Tools::GFF;
use BILS::Handler::GXFhandler qw(:Ok);
use BILS::GFF3::Statistics qw(:Ok);
use BILS::Handler::GFF3handler qw(:Ok);
# END libraries

my $header = qq{
########################################################
# BILS 2015 - Sweden                                   #  
# jacques.dainat\@bils.se                               #
# Please cite BILS (www.bils.se) when using this tool. #
########################################################
};

# PARAMETERS - OPTION
my $opt_reffile;
my $opt_output;
my $opt_BlastFile;
my $opt_InterproFile;
my $opt_name;
my $opt_nameU;
my $optFillFrame;
my $optForceFillFrame;
my $optOnlyStat;
my $opt_genomeSize;
my $opt_removeUTR;
my $opt_removemRNAduplicated;
my $opt_help = 0;
# END PARAMETERS - OPTION

# for ID name
my $nbGeneName;
my $nbmRNAname;
my $nbCDSname;
my $nbExonName;
my $nbUTRName;
my $nbRepeatName;
# END ID name

# FOR FUNCTIONS BLAST#
my %nameBlast;
my %geneNameBlast;
my %mRNANameBlast;
my %mRNAproduct;
my %geneNameGiven;
my %duplicateNameGiven;
my $nbDuplicateNameGiven=0;
my $nbDuplicateName=0;
my $nbNamedGene=0;
my $nbGeneNameInBlast=0;
# END FOR FUNCTION BLAST#

# FOR FUNCTIONS INTERPRO#
my %TotalTerm;
my %GeneAssociatedToTerm;
my %mRNAAssociatedToTerm;
my %functionData;
my %functionStreamOutput;
my %geneWithoutFunction;
my %geneWithFunction;
my $nbmRNAwithoutFunction=0;
my $nbmRNAwithFunction=0;
my $nbGeneWithGOterm=0;
my $nbTotalGOterm=0;
# END FOR FUNCTION INTERPRO#

# OPTION MANAGMENT
my @copyARGV=@ARGV;
if ( !GetOptions( 'f|ref|reffile|gff|gff3=s' => \$opt_reffile,
                  'b|blast=s' => \$opt_BlastFile,
                  'i|interpro=s' => \$opt_InterproFile,
                  'id=s' => \$opt_name,
                  'idau=s' => \$opt_nameU,               
                  'g|gs=s' => \$opt_genomeSize,
                  'gf=i'      => \$nbGeneName,
                  'mf=i'      => \$nbmRNAname,
                  'cf=i'      => \$nbCDSname,
                  'ef=i'      => \$nbExonName,
                  'uf=i'      => \$nbUTRName,
                  'rf=i'      => \$nbRepeatName,
                  'ff'      => \$optFillFrame,
                  's'      => \$optOnlyStat,
                  'o|output=s'      => \$opt_output,

                  'h|help!'         => \$opt_help ) )
{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

if ($opt_help) {
    pod2usage( { -verbose => 2,
                 -exitval => 2,
                 -message => "$header\n" } );
}

if ( ! (defined($opt_reffile)) ){
    pod2usage( {
           -message => "$header\nAt least 1 parameter is mandatory:\nInput reference gff file (--f)\n\n".
           "Many optional parameters are available. Look at the help documentation to know more.\n",
           -verbose => 0,
           -exitval => 1 } );
}

# counters for ids initialisation
if (! $nbGeneName){$nbGeneName=1};
if (! $nbmRNAname){$nbmRNAname=1};
if (! $nbCDSname){$nbCDSname=1};
if (! $nbExonName){$nbExonName=1};
if (! $nbUTRName){$nbUTRName=1};
if (! $nbRepeatName){$nbRepeatName=1};


#################################################
####### START Manage files (input output) #######
#################################################

my $streamBlast = IO::File->new();
my $streamInter = IO::File->new();

# Manage Blast File
if (defined $opt_BlastFile){
$streamBlast->open( $opt_BlastFile, 'r' ) or
  croak(
     sprintf( "Can not open '%s' for reading: %s", $opt_BlastFile, $! ) );
}

# Manage Interpro file
if (defined $opt_InterproFile){
$streamInter->open( $opt_InterproFile, 'r' ) or
  croak(
     sprintf( "Can not open '%s' for reading: %s", $opt_InterproFile, $! ) );
}

##########################
##### Manage Output ######
my @outputTab;

if (defined($opt_output) ) {
  if (-f $opt_output){
      print "Cannot create a directory with the name $opt_output because a file with this name already exists.\n";exit();
  }
  if (-d $opt_output){
      print "The output directory choosen already exists. Please geve me another Name.\n";exit();
  }
  #### Case 1 => option ouput option onlyStat
  mkdir $opt_output;

  my $ostreamReport=IO::File->new(">".$opt_output."/report.txt" ) or
  croak( sprintf( "Can not open '%s' for writing %s", $opt_output."/duplicatedIDFeatures.gff", $! ));
  push (@outputTab, $ostreamReport);

  #### Case 2 => option ouput NO option onlyStat
  if (! $optOnlyStat){
      my $ostreamCoding=Bio::Tools::GFF->new(-file => ">".$opt_output."/AllFeatures.gff", -gff_version => 3 ) or
      croak(sprintf( "Can not open '%s' for writing %s", $opt_output."AllFeatures.gff", $! ));
      push (@outputTab, $ostreamCoding);
      
      my $ostreamNormalGene=Bio::Tools::GFF->new(-file => ">".$opt_output."/codingGeneFeatures.gff", -gff_version => 3 ) or
      croak( sprintf( "Can not open '%s' for writing %s", $opt_output."/codingGeneFeatures.gff", $! ));
      push (@outputTab, $ostreamNormalGene);

      my $ostreamOtherRNAGene=Bio::Tools::GFF->new(-file => ">".$opt_output."/otherRNAfeatures.gff", -gff_version => 3 ) or
      croak(sprintf( "Can not open '%s' for writing %s", $opt_output."/otherRNAfeatures.gff", $! ));
      push (@outputTab, $ostreamOtherRNAGene);

      my $ostreamRepeats=Bio::Tools::GFF->new(-file => ">".$opt_output."/repeatsFeatures.gff", -gff_version => 3 )or
      croak( sprintf( "Can not open '%s' for writing %s", $opt_output."/repeatsFeatures.gff", $! ));
      push (@outputTab, $ostreamRepeats);
  }
}
### Case 3 => No output option => everithing will be display on screen. 
### Case 4 => If option onlyStat provided the script will stop before writting results.
else {
 my $ostreamReport = \*STDOUT or die ( sprintf( "Can not open '%s' for writing %s", "STDOUT", $! ));
  push (@outputTab, $ostreamReport);

  my $ostream  = IO::File->new();
  $ostream->fdopen( fileno(STDOUT), 'w' ) or croak( sprintf( "Can not open STDOUT for writing: %s", $! ) );
  my $outputGFF = Bio::Tools::GFF->new( -fh => $ostream ) or croak( sprintf( "Can not open STDOUT for writing: %s", $! ) );

  #my $outputGFF = Bio::Tools::GFF->new( \*STDOUT, -gff_version => 3 ) or
  #croak( sprintf( "Can not open STDOUT for writing: %s", $! ) );
  push (@outputTab, $outputGFF);
  push (@outputTab, $outputGFF);
  push (@outputTab, $outputGFF);
  push (@outputTab, $outputGFF);
  push (@outputTab, $outputGFF); ### Creation of a list of output stream <= In this case every time the same ! Because it for display to the screen                                 
}

###############################################
####### END Manage files (input output) #######
###############################################
my $stringPrint = strftime "%m/%d/%Y at %Hh%Mm%Ss", localtime;

$stringPrint .= "\nusage: $0 @copyARGV\n";
$stringPrint .= "vvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n";
$stringPrint .= "vvvvvvvv OPTION INFO vvvvvvvv\n\n";
$stringPrint .= "->We will calculate statistics about the input file\n";
if ($optOnlyStat){ $stringPrint .= "->We will just calculate the statistics. Features will not be printed.\n";}
else{ $stringPrint .= "->The feature will be sorted before to print them\n";}

my $prefixName;
if ($opt_name){
  $prefixName=$opt_name;
  $stringPrint .= "->IDs will be changed using $opt_name as prefix.\nIn the case of discontinuous features (i.e. a single feature that exists over multiple genomic locations) the same ID may appear on multiple lines.".
  " All lines that share an ID collectively represent a signle feature.\n";
  $stringPrint .= "-> Exon will be expanded even if not asked to avoid loss of multiple parent during exon renaming.\n"
}
if ($opt_nameU){
  $stringPrint .= "->IDs will be changed using $opt_nameU as prefix. Id of features that share an ID collectively will be change in different and uniq ID.\n";
  $prefixName=$opt_nameU;
}
if($optFillFrame or $optForceFillFrame){
  $stringPrint .= "->CDS frame will be fill\n";
  $stringPrint .= "-> Exon will be expanded even if not asked to avoid loss of multiple parent during exon renaming.\n"
}
$stringPrint .= "\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n\n";

# Display
$outputTab[0]->print($stringPrint);
if($opt_output){print "$stringPrint";} # When ostreamReport is a file we have to also display on screen



                                                      #######################
                                                      #        MAIN         #
#                     >>>>>>>>>>>>>>>>>>>>>>>>>       #######################       <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) = BILS::Handler::GFF3handler->slurp_gff3_file_JD($opt_reffile);
print("Parsing Finished\n\n");
### END Parse GFF input #
#########################

#Print directly what has been read 
#my $stat = gff3_statistics($hash_omniscient, $opt_genomeSize);
#  foreach my $info (@$stat){
#    $outputTab[0]->print("$info");
#  }

################################
# MANAGE FUNCTIONAL INPUT FILE #

# Manage Blast File
if (defined $opt_BlastFile){
  parseAnnieBlast($streamBlast,$opt_BlastFile);
}

# Manage Interpro File
if (defined $opt_InterproFile){
  parseAnnieInterpro($streamInter,$opt_InterproFile);
  
  # create streamOutput
  if($opt_output){
    foreach my $type (keys %functionData){
      my $ostreamFunct = IO::File->new(); 
      $ostreamFunct->open( $opt_output."/$type.txt", 'w' ) or
          croak(
              sprintf( "Can not open '%s' for writing %s", $opt_output."/$type.txt", $! )
          );
      $functionStreamOutput{$type}=$ostreamFunct;
    }
  }
}
# END MANAGE FUNCTIONAL INPUT FILE #
####################################

#################################
# GO THROUGH OMISCIENT          # Will create

my %omniscient_gene;
my %omniscient_other;
my %omniscient_repeat;
my @list_geneID_l1;
my @list_OtherRnaID_l1;
my @list_repeatID_l1;
#################
# create list by 3 type of feature (gene, trna, repeats). Allows to create different outputs

# level 1
foreach my $primary_tag_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_level1 = gene or repeat etc...
  foreach my $id_level1 (keys %{$hash_omniscient->{'level1'}{$primary_tag_level1}}){
    if($primary_tag_level1 =~ /repeat/){
      push(@list_repeatID_l1, $id_level1)
    }
    else{
      # get one level2 feature to check wich level1 feature it is
      foreach my $primary_tag_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_level2 = mrna or mirna or ncrna or trna etc...        
        if ( exists ($hash_omniscient->{'level2'}{$primary_tag_level2}{$id_level1} ) ){
          my $one_feat=@{$hash_omniscient->{'level2'}{$primary_tag_level2}{$id_level1}}[0];
          if(lc($one_feat->primary_tag) eq "mrna"){
            push(@list_geneID_l1, $id_level1);
            last;
          }
          else{
            push(@list_OtherRnaID_l1, $id_level1);
            last;
          }
        }
      }
    }
  }
}

##########################
# create sub omniscients
my %hash_of_omniscient;
if(@list_geneID_l1){
  fill_omniscient_from_other_omniscient_level1_id(\@list_geneID_l1, $hash_omniscient, \%omniscient_gene);
  $hash_of_omniscient{'Coding_Gene'}=\%omniscient_gene;
}
if(@list_OtherRnaID_l1){
  fill_omniscient_from_other_omniscient_level1_id(\@list_OtherRnaID_l1, $hash_omniscient, \%omniscient_other);
  $hash_of_omniscient{'Non_Coding_Gene'}=\%omniscient_other;
}
if(@list_repeatID_l1){
  fill_omniscient_from_other_omniscient_level1_id(\@list_repeatID_l1, $hash_omniscient, \%omniscient_repeat);
  $hash_of_omniscient{'Repeat'}=\%omniscient_repeat;
}

##############
# STATISTICS #
foreach my $key_hash (keys %hash_of_omniscient){
  $outputTab[0]->print("Information about $key_hash\n");
  if($opt_output){print "Information about $key_hash\n";} # When ostreamReport is a file we have to also display on screen
  my $hash_ref = $hash_of_omniscient{$key_hash};
  my $stat;
  if($opt_genomeSize){
    $stat = gff3_statistics($hash_ref, $opt_genomeSize);
  }else{$stat = gff3_statistics($hash_ref);}

  #print statistics
  foreach my $infoList (@$stat){
    foreach my $info (@$infoList){
      $outputTab[0]->print("$info");
      if($opt_output){print "$info";} # When ostreamReport is a file we have to also display on screen
    }
    $outputTab[0]->print("\n");
    if($opt_output){print "\n";} # When ostreamReport is a file we have to also display on screen
  }
}

# END STATISTICS #
##################

# STOP script here if option given
if ($optOnlyStat) {
    print "only statistics option => We stop here.\nEND";exit;
}

###################
#Fil frame is asked
foreach my $key_hash (keys %hash_of_omniscient){
  my $hash_ref = $hash_of_omniscient{$key_hash};
  if($optFillFrame){
    fil_cds_frame($key_hash);
  }
}

###########################
# change FUNCTIONAL information if asked for
if ($opt_BlastFile || $opt_InterproFile ){#|| $opt_BlastFile || $opt_InterproFile){
    my $hash_ref = $hash_of_omniscient{'Coding_Gene'};

    #################
    # == LEVEL 1 == #
    #################
    foreach my $primary_tag_level1 (keys %{$hash_ref ->{'level1'}}){ # primary_tag_level1 = gene or repeat etc...
      foreach my $id_level1 (keys %{$hash_ref ->{'level1'}{$primary_tag_level1}}){
        
        my $feature_level1=$hash_ref->{'level1'}{$primary_tag_level1}{$id_level1};
        # Clean NAME attribute
        $feature_level1->remove_tag('Name');

        #Manage Name if otpion setting
        if( $opt_BlastFile ){
          if (exists ($geneNameBlast{$id_level1})){
            create_or_replace_tag($feature_level1, 'Name', $geneNameBlast{$id_level1});
            $nbNamedGene++;
            
            # Check name duplicated given
            my $nameClean=$geneNameBlast{$id_level1};
            $nameClean =~ s/_([2-9]{1}[0-9]*|[0-9]{2,})*$//;
            
            my $nameToCompare;
            if(exists ($nameBlast{$nameClean})){ # We check that is really a name where we added the suffix _1
              $nameToCompare=$nameClean;
            }
            else{$nameToCompare=$geneNameBlast{$id_level1};} # it was already a gene_name like BLABLA_12

            if(exists ($geneNameGiven{$nameToCompare})){
                $nbDuplicateNameGiven++; # track total
                $duplicateNameGiven{$nameToCompare}++; # track diversity
            }
            else{$geneNameGiven{$nameToCompare}++;} # first time we have given this name
          }
        }

        #################
        # == LEVEL 2 == #
        #################
        foreach my $primary_tag_key_level2 (keys %{$hash_ref->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
          
          if ( exists ($hash_ref->{'level2'}{$primary_tag_key_level2}{$id_level1} ) ){
            foreach my $feature_level2 ( @{$hash_ref->{'level2'}{$primary_tag_key_level2}{$id_level1}}) {

              my $level2_ID = lc($feature_level2->_tag_value('ID'));
              # Clean NAME attribute
              $feature_level2->remove_tag('Name');

              #Manage Name if option set
              if($opt_BlastFile){
                if (exists ($mRNANameBlast{$level2_ID})){
                  my $mRNABlastName=$mRNANameBlast{$level2_ID};
                  create_or_replace_tag($feature_level2, 'Name', $mRNABlastName);
                }
                my $productData=printProductFunct($level2_ID);
                if ($productData ne ""){
                  create_or_replace_tag($feature_level2, 'description', $productData);
                }
                else {
                  create_or_replace_tag($feature_level2, 'description', "hypothetical protein");
                } #Case where the protein is not known
              }

              # print function if option
              if($opt_InterproFile){
                my $parentID=$feature_level2->_tag_value('Parent');

                if (addFunctions($feature_level2, $opt_output)){
                  $nbmRNAwithFunction++;$geneWithFunction{$parentID}++;
                  if(exists ($geneWithoutFunction{$parentID})){
                    delete $geneWithoutFunction{$parentID};
                  }
                }
                else{
                  $nbmRNAwithoutFunction++;
                  if(! exists ($geneWithFunction{$parentID})){
                    $geneWithoutFunction{$parentID}++;
                  }
                }
              }
            }
          }
        }
      } 
    }
}


###########################
# change names if asked for
if ($opt_nameU || $opt_name ){#|| $opt_BlastFile || $opt_InterproFile){
  foreach my $key_hash (keys %hash_of_omniscient){
    my $hash_ref = $hash_of_omniscient{$key_hash};

    #################
    # == LEVEL 1 == #
    #################
    foreach my $primary_tag_level1 (keys %{$hash_ref ->{'level1'}}){ # primary_tag_level1 = gene or repeat etc...
      foreach my $id_level1 (keys %{$hash_ref ->{'level1'}{$primary_tag_level1}}){
        
        my $feature_level1=$hash_ref->{'level1'}{$primary_tag_level1}{$id_level1};
        my $newID_level1;

        #keep track of Maker ID
        if($opt_BlastFile){#In that case the name given by Maker is removed from ID and from Name. We have to kee a track
          create_or_replace_tag($feature_level1, 'makerName', $id_level1);
        }

        if(lc($primary_tag_level1) =~ /repeat/ ){
          $newID_level1 = manageID($prefixName,$nbRepeatName,'R'); 
          $nbRepeatName++;
          create_or_replace_tag($feature_level1, 'ID', $newID_level1);
        }
        else{
          $newID_level1 = manageID($prefixName,$nbGeneName,'G'); 
          $nbGeneName++; 
          create_or_replace_tag($feature_level1, 'ID', $newID_level1);
        }

        #################
        # == LEVEL 2 == #
        #################
        foreach my $primary_tag_key_level2 (keys %{$hash_ref->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
          
          if ( exists ($hash_ref->{'level2'}{$primary_tag_key_level2}{$id_level1} ) ){
            foreach my $feature_level2 ( @{$hash_ref->{'level2'}{$primary_tag_key_level2}{$id_level1}}) {

              my $level2_ID = lc($feature_level2->_tag_value('ID'));
              my $newID_level2;
              
              #keep track of Maker ID
              if($opt_InterproFile){#In that case the name given by Maker is removed from ID and from Name. We have to kee a track
                create_or_replace_tag($feature_level2, 'makerName', $level2_ID);
              }

              if(lc($feature_level2) =~ /repeat/ ){
                print "What should we do ? implement something. L1 and l2 repeats will have same name ...\n";exit;
              }
              else{
                $newID_level2 = manageID($prefixName,$nbmRNAname,"T");
                $nbmRNAname++; 
                create_or_replace_tag($feature_level2, 'ID', $newID_level2);
                create_or_replace_tag($feature_level2, 'Parent', $newID_level1);
              }
        
              #################
              # == LEVEL 3 == #
              #################
             
              foreach my $primary_tag_level3 (keys %{$hash_ref->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
                
                  if ( exists ($hash_ref->{'level3'}{$primary_tag_level3}{$level2_ID} ) ){

                    foreach my $feature_level3 ( @{$hash_ref->{'level3'}{$primary_tag_level3}{$level2_ID}}) {
                      
                      #keep track of Maker ID
                      my $level3_ID = lc($feature_level2->_tag_value('ID'));
                      if($opt_InterproFile){#In that case the name given by Maker is removed from ID and from Name. We have to kee a track
                        create_or_replace_tag($feature_level3, 'makerName', $level3_ID);
                      }

                      if($primary_tag_level3 =~ /exon/ ){
                        my $newID = manageID($prefixName,$nbExonName,'E'); 
                          $nbExonName++;
                          create_or_replace_tag($feature_level3, 'ID', $newID);
                          create_or_replace_tag($feature_level3, 'Parent', $newID_level2);

                      }
                      elsif($primary_tag_level3 =~ /cds/){
                        my $newID = manageID($prefixName,$nbCDSname,'C'); 
                        if($opt_nameU){$nbCDSname++;}
                        create_or_replace_tag($feature_level3, 'ID', $newID);
                        create_or_replace_tag($feature_level3, 'Parent', $newID_level2);
                      }

                      elsif($primary_tag_level3 =~ /utr/){
                        my $newID = manageID($prefixName,$nbUTRName,'U');
                        if($opt_nameU){$nbUTRName++;}
                        create_or_replace_tag($feature_level3, 'ID', $newID);
                        create_or_replace_tag($feature_level3, 'Parent', $newID_level2);
                      }
                      else{
                        print "What that one ? $primary_tag_level3\n";
                      }
                      push (@{$hash_ref->{'level3'}{$primary_tag_level3}{lc($newID_level2)}}, $feature_level3);
                    }
                  }
                  if ($opt_name and  $primary_tag_level3 =~ /utr/){$nbUTRName++;} # with this option we increment UTR name only for each UTR 
                  if ($opt_name and  $primary_tag_level3 =~ /cds/){$nbCDSname++;} # with this option we increment cds name only for each cds 
              }
            }
          }
        }
      }
    }
  }
}

###########################
# RESULT PRINTING
###########################



##############################
# print FUNCITONAL INFORMATION
$stringPrint =""; # reinitialise (use at the beginning)
if ($opt_InterproFile){
  #print INFO
  my $lineB=       "___________________________________________________________________________________________________";
  $stringPrint .= " ".$lineB."\n";
  $stringPrint .= "|          | Nb Total term | Nb mRNA with term  | Nb mRNA updated by term | Nb gene updated by term |\n";
  $stringPrint .= "|          | in Annie File |   in Annie File    | in our annotation file  | in our annotation file  |\n";
  $stringPrint .= "|".$lineB."|\n";

  foreach my $type (keys %functionData){
    my $total_type = $TotalTerm{$type};
    my $mRNA_type_Annie = keys %{$functionData{$type}};
    my $mRNA_type = keys %{$mRNAAssociatedToTerm{$type}};
    my $gene_type = keys %{$GeneAssociatedToTerm{$type}};
    $stringPrint .= "|".sizedPrint(" $type",10)."|".sizedPrint($total_type,15)."|".sizedPrint($mRNA_type_Annie,20)."|".sizedPrint($mRNA_type,25)."|".sizedPrint($gene_type,25)."|\n|".$lineB."|\n";
  }

  #RESUME TOTAL OF FUNCTION ATTACHED
  my $listOfFunction;
  foreach my $funct (keys %functionData){
    $listOfFunction.="$funct,";
  }
  chop $listOfFunction;
  $stringPrint .= "nb mRNA without Functional annotation ($listOfFunction) = $nbmRNAwithoutFunction\n";
  $stringPrint .= "nb mRNA with Functional annotation ($listOfFunction) = $nbmRNAwithFunction\n";
  my $nbGeneWithoutFunction= keys %geneWithoutFunction;
  $stringPrint .= "nb gene without Functional annotation ($listOfFunction) = $nbGeneWithoutFunction\n";
  my $nbGeneWithFunction= keys %geneWithFunction;
  $stringPrint .= "nb gene with Functional annotation ($listOfFunction) = $nbGeneWithFunction\n";
  
}

if($opt_BlastFile){
  my $nbGeneDuplicated=keys %duplicateNameGiven;
  $nbDuplicateNameGiven=$nbDuplicateNameGiven+$nbGeneDuplicated; # Until now we have counted only name in more, now we add the original name.
  $stringPrint .= "$nbGeneNameInBlast gene names have been retrieved in the blast file. Among them there are $nbDuplicateName gene names duplicated. \n".
  "$nbNamedGene genes have been named. $nbGeneDuplicated names are shared at least per two genes for a total of $nbDuplicateNameGiven genes.\n";

  if($opt_output){
    my $duplicatedNameOut=IO::File->new(">".$opt_output."/duplicatedNameFromBlast.txt" );
    foreach my $name (sort { $duplicateNameGiven{$b} <=> $duplicateNameGiven{$a} } keys %duplicateNameGiven){
      print $duplicatedNameOut "$name\t".($duplicateNameGiven{$name}+1)."\n";
    }
  }
}


# Display
$outputTab[0]->print("$stringPrint");
if(defined $opt_output){print "$stringPrint";}

####################
# PRINT IN FILES
####################
#print step
printf("Writing result\n");
#print gene (mRNA)
print_omniscient(\%omniscient_gene, $outputTab[1]);
print_omniscient(\%omniscient_gene, $outputTab[2]);
#print other RNA gene
print_omniscient(\%omniscient_other, $outputTab[1]);
print_omniscient(\%omniscient_other, $outputTab[3]);
#print repeat
print_omniscient(\%omniscient_repeat, $outputTab[1]);
print_omniscient(\%omniscient_repeat, $outputTab[4]);

      ######################### 
      ######### END ###########
      #########################
#######################################################################################################################
        ####################
         #     methods    #
          ################
           ##############
            ############
             ##########
              ########
               ######
                ####
                 ##

# each mRNA of a gene has its proper gene name. Most often is the same, and annie added a number at the end. To provide only one gene name, we remove this number and then remove duplicate name (case insensitive).
# If it stay at the end of the process more than one name, they will be concatenated together.
# It remove redundancy intra name.
sub manageGeneNameBlast{
  my ($geneName)=@_;
  foreach my $element (keys %$geneName){
    my @tab=@{$geneName->{$element}};
    
    my %seen;
    my @unique;
    for my $w (@tab) { # remove duplicate in list case insensitive
      $w =~ s/_[0-9]+$// ;
      next if $seen{lc($w)}++;
      push(@unique, $w);
    }

    my $finalName="";
    my $cpt=0;
    foreach my $name (@unique){  #if several name we will concatenate them together
        if ($cpt == 0){
          $finalName .="$name";
        }
        else{$finalName .="_$name"}
    }
    $geneName->{$element}=$finalName;
    $nameBlast{lc($finalName)}++;
  }
}

# creates gene ID correctly formated (PREFIX,TYPE,NUMBER) like HOMSAPG00000000001 for a Homo sapiens gene.
sub manageID{
  my ($prefix,$nbName,$type)=@_;
  my $result="";
  my $numberNum=11;
  my $GoodNum="";
  for (my $i=0; $i<$numberNum-length($nbName); $i++){
    $GoodNum.="0";
  }
  $GoodNum.=$nbName;
  $result="$prefix$type$GoodNum";

  return $result;
}

# Create String containing the product information associated to the mRNA
sub printProductFunct{
  my ($refname)=@_;
  my $String="";
  my $first="yes";
  if (exists $mRNAproduct{$refname}){
    foreach my $element (@{$mRNAproduct{$refname}})
    { 
      if($first eq "yes"){
        $String.="$element";
        $first="no";
      }
      else{$String.=",$element";}
    }
  }
  return $String;
}

sub addFunctions{
  my ($feature, $opt_output)=@_;

  my $functionAdded=undef;
  my $ID=lc($feature->_tag_value('ID'));
  foreach my $function_type (keys %functionData){
    
    
    if(exists ($functionData{$function_type}{$ID})){
      $functionAdded="true";

      my $data_list;

      if(lc($function_type) eq "go"){
        foreach my $data (@{$functionData{$function_type}{$ID}}){
          $feature->add_tag_value('Ontology_term', $data);
          $data_list.="$data,";
        }
      }
      else{
        foreach my $data (@{$functionData{$function_type}{$ID}}){
          $feature->add_tag_value('Dbxref', $data);
          $data_list.="$data,";
        }
      }

      if ($opt_output){
          my $streamOut=$functionStreamOutput{$function_type};
          my $ID = $feature->_tag_value('ID');
          chop $data_list;
          print $streamOut "$ID\t$data_list\n";
        }
    }
  }
  return $functionAdded;
}

# method to par annie blast file
sub parseAnnieBlast {
  my($file_in,$fileName) = @_;
  print( "Reading features from $fileName...\n");
  my %geneName; my %linkBmRNAandGene;
  my $geneID ="";
  my $nameGene="";

  #FIRST PARSE THE FILE
  while( my $line = <$file_in>)  {   
    next if $line =~ /^\s*$/;
    my @values = split(/\s/, $line);

    if ($values[1] eq "name"){
      $geneID=$values[0];
      $nameGene=$values[2];
      $nameGene=~ s/\n//g;
      push ( @{ $geneName{lc($geneID)} }, lc($nameGene) );
    }
    elsif ($values[1] eq "product"){
      my $mRNAID=$values[0];
      @values = split(/\sproduct\s/, $line);
      my $product=$values[1];
      $product=~ s/\n//g;
      push ( @{ $mRNAproduct{lc($mRNAID)} }, $product );
      if ($nameGene ne ""){
        push( @{ $linkBmRNAandGene{lc($geneID)}}, lc($mRNAID)); # save mRNA name for each gene name 
        $geneID ="";
        $nameGene="";
      }
    }
    else {print "/!\\ Achtung !! something strange in this file ... line is: $line";}   
  }

  # secondly Manage NAME (If several)
  manageGeneNameBlast(\%geneName); # Remove redundancy to have only one name for each gene
  
  #Then CLEAN NAMES REDUNDANCY inter gene
  my %geneNewNameUsed;
  foreach my $geneID (keys %geneName){

    $nbGeneNameInBlast++;
    
    my @mRNAList=@{$linkBmRNAandGene{$geneID}};
    my $String = $geneName{$geneID};
#    print "$String\n";
    if (! exists( $geneNewNameUsed{$String})){
      $geneNewNameUsed{$String}++;
      $geneNameBlast{$geneID}=$String;
      # link name to mRNA and and isoform name _1 _2 _3 if several mRNA
      my $cptmRNA=1;
      if ($#mRNAList != 0) {
        foreach my $mRNA (@mRNAList){
          $mRNANameBlast{$mRNA}=$String."_iso".$cptmRNA;
          $cptmRNA++;
        }
      }
      else{$mRNANameBlast{$mRNAList[0]}=$String;}
    }
    else{ #in case where name was already used, we will modified it by addind a number like "_2"
      $nbDuplicateName++;
      $geneNewNameUsed{$String}++;
      my $nbFound=$geneNewNameUsed{$String};
      $String.="_$nbFound";
      $geneNewNameUsed{$String}++;
      $geneNameBlast{$geneID}=$String;
      # link name to mRNA and and isoform name _1 _2 _3 if several mRNA
      my $cptmRNA=1;  
      if ($#mRNAList != 0) {
        foreach my $mRNA (@mRNAList){
          $mRNANameBlast{$mRNA}=$String."_iso".$cptmRNA;
          $cptmRNA++;
        }
      }
      else{$mRNANameBlast{$mRNAList[0]}=$String;}
    }
  }
}

# method to par annie Interpro file
sub parseAnnieInterpro {
  my($file_in,$fileName) = @_;
  print( "Reading features from $fileName...\n");
  
  while( my $line = <$file_in>)  {    
    my @values = split(/\t/, $line);
    my $mRNAID=lc($values[0]);
    if ((lc($values[1]) eq "db_xref") || (lc($values[1]) eq "dbxref")){
      my $data=$values[2];
      $data=~ s/\n//g;
      $data=~ s/\s//g; #remove space
      my @element = split(/:/,$data);
      my $typeEl = $element[0];
      @element = split(/\|/,$data); #cut at character | 

      foreach my $oneEl (@element){
        $TotalTerm{$typeEl}++;
        push ( @{$functionData{$typeEl}{$mRNAID}} , $oneEl );
        if ( exists $hash_mRNAGeneLink->{$mRNAID}){ ## check if exists among our current gff annotation file analyzed
          $mRNAAssociatedToTerm{$typeEl}{$mRNAID}++;
          $GeneAssociatedToTerm{$typeEl}{$hash_mRNAGeneLink->{$mRNAID}}++;
        }
      }        
    }      
  }
}

sub sizedPrint{
  my ($term,$size) = @_;
  my $result; my $sizeTerm=length($term);
  if ($sizeTerm > $size ){
    $result=substr($term, 0,$size);
    return $result;
  }
  else{
    my $nbBlanc=$size-$sizeTerm;
    $result=$term;
    for (my $i = 0; $i < $nbBlanc; $i++){
      $result.=" ";
    }
    return $result;
  }
}

__END__

=head1 NAME

gff3manager_JD.pl -
The script take a gff3 file as input. -
Without option the script only sort the data. -
With corresponding parameters, it can add functional annotations from <annie> output files
>The blast against Prot Database file from annie allows to fill the field NAME for gene and PRODUCT for mRNA.
>The blast against Interpro Database tsv file from annie allows to fill the DBXREF field with pfam, tigr, interpro and GO terms data.
The script expand exons sharing multiple mRNA (Parent attributes contains multiple parental mRNA). One exon by parental mRNA will be created.
With the <id> option the script will change all the ID field by an Uniq ID created from the given prefix, a letter to specify the kind of feature (G,T,C,E,U), and the feature number.

The result is written to the specified output file, or to STDOUT.
Remark: If there is duplicate in the file they will be removed in the output. In that case you should be informed.

=head1 SYNOPSIS

    ./gff3manager_JD.pl -f=infile.gff [ -b annie_blast_infile -i annie_interpro_infile.tsv -x annie_interpro_infile.xml -e --id ABCDEF [-gf 20] -s -utr -utrr 10 --output outfile ]
    ./gff3manager_JD.pl --help

=head1 OPTIONS

=over 8

=item B<-f>, B<--reffile>,B<-ref> , B<--gff> or B<--gff3> 

Input GFF3 file that will be read (and sorted)

=item B<-b> or B<--blast> 

Input annie blast file that will be used to complement the features read from
the first file (specified with B<--ref>).

=item B<-i> or B<--interpro> 

Input annie interpro file (.tsv) that will be used to complement the features read from
the first file (specified with B<--ref>).

=item B<-g> or B<--gs> 

This option inform about the genome size in oder to compute more statistics. You can give the size in Nucleotide or directly the fasta file.

=item B<-id>

This option will changed the id name. It will create from id prefix (usually 6 letters) given as input, uniq IDs like prefixE00000000001. Where E mean exon. Instead E we can have C for CDS, G for gene, T for mRNA, U for Utr.
In the case of discontinuous features (i.e. a single feature that exists over multiple genomic locations) the same ID may appear on multiple lines. All lines that share an ID collectively represent a signle feature.

=item B<-idau>

This option (id all uniq) is similar to -id option but Id of features that share an ID collectively will be change by different and uniq ID.

=item B<-gf>

Usefull only if -id is used.
This option is used to define the number that will be used to begin to number the gene id (gf for "gene from"). By default begin by 1.

=item B<-mf>

Usefull only if -id is used.
This option is used to define the number that will be used to begin to number the mRNA id (mf for "mRNA from"). By default begin by 1.

=item B<-cf>

Useful only if -id is used.
This option is used to define the number that will be used to begin to number the CDS id (cf for "CDS from"). By default begin by 1.

=item B<-ef>

Useful only if -id is used.
This option is used to define the number that will be used to begin to number the exon id (ef for "Exon from"). By default begin by 1.

=item B<-uf>

Useful only if -id is used.
This option is used to define the number that will be used to begin to number the UTR id (uf for "UTR from"). By default begin by 1.

=item B<-rf>

Useful only if -id is used.
This option is used to define the number that will be used to begin to number the repeat id (rf for "Repeat from"). By default begin by 1.

=item B<-ff>

ff means fill frame.
This option is used to add the CDS frame. If frames already exist, the script overwrite them.

=item B<-s>

Just compute some statisctics about the input file and stop.

=item B<-o> or B<--output>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.

=back

=cut
