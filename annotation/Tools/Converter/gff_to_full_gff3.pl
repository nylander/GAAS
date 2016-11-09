#!/usr/bin/env perl

## TO DO:
## Dont manage if some isoform (not all) should be modelate for a same gene. (gene will appear 2 times)
# If tehre is only the level2 to build

use strict;
use Getopt::Long;
use Clone 'clone';
use Pod::Usage;
use Bio::Tools::GFF;
use BILS::Handler::GTFhandler qw(:Ok);
use BILS::Handler::GFF3handler qw(:Ok);
use BILS::Handler::GXFhandler qw(:Ok);
use Data::Dumper;

my $header = qq{
########################################################
# BILS 2015 - Sweden                                   #  
# jacques.dainat\@bils.se                               #
# Please cite BILS (www.bils.se) when using this tool. #
########################################################
};

my $cptgg;
my $outfile = undef;
my $gff = undef;
my $attributes=undef;
my $keepAllAtt=undef;
my $help= 0;

if ( !GetOptions(
    "help|h" => \$help,
    "gff|in|infile=s" => \$gff,
    "attributes|attribute|a|att=s" => \$attributes,
    "kaa" => \$keepAllAtt,
    "outfile|out|o|output=s" => \$outfile))

{
    pod2usage( { -message => "Failed to parse command line.",
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($help) {
    pod2usage( { -message => "$header", 
                 -verbose => 2,
                 -exitval => 2 } );
}
 
if ( ! (defined($gff)) ){
    pod2usage( {
           -message => "\nAt least 1 parameter is mandatory:\nInput reference gff file\n",
           -verbose => 0,
           -exitval => 1 } );
}

## Manage output file
my $gffout;
if ($outfile) {
open(my $fh, '>', $outfile) or die "Could not open file '$outfile' $!";
  $gffout= Bio::Tools::GFF->new(-fh => $fh, -gff_version => 3 );
}
else{
  $gffout = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 3);
}

# Mange attributes if given
### If attributes given, parse them:
my %attListOk;
my @attListPair;
if ($attributes){
  if($attributes != 'all'){
    @attListPair= split(/,/, $attributes);
    my $nbAtt=$#attListPair+1;
    print "In addition of gene_id attribute, we will keep $nbAtt attribute(s):\n";
    foreach my $attributeTuple (@attListPair){ 
      my @attList= split(/\//, $attributeTuple);
      if($#attList == 0){ # Attribute alone
        $attListOk{$attList[0]}='null';
        print "$attList[0]\n";
      }
      else{ # Attribute we have to replace by a new name
        $attListOk{$attList[0]}=$attList[1];
        print "$attList[0] (Will be replaced by $attList[1])\n";
      }
    }
    print "\n";
  }
  else{print "You chose to keep all attributes !\n";}
}



### Parse GFF input file and add annotations
my ($hash_omniscient, $hash_mRNAGeneLink) = BILS::Handler::GXFhandler->slurp_gff3_file_JD($gff);

## Check duplication
#my $size = keys %$hash_duplicatedFeature;
#if($size == 0){
#	print "No duplicated feature found !\n";
#}else { print "//// !! ACHTUNG \\\\\\\\\ $size duplicated feature found !";}

my $hash_level1=$hash_omniscient->{"level1"};
my $hash_level2=$hash_omniscient->{"level2"};
my $hash_level3=$hash_omniscient->{"level3"};
my %Level2featuresNames;
my %hash_sortBySeq;

### check omniscient
my @keysLevel1 = keys %{$hash_level1};
my $sizeLevel1 = @keysLevel1;
my @keysLevel2 = keys %{$hash_level2};
my $sizeLevel2 = @keysLevel2;
my @keysLevel3 = keys %{$hash_level3};
my $sizeLevel3 = @keysLevel3;

my $sizeLevels = $sizeLevel1+$sizeLevel2+$sizeLevel3;
if($sizeLevels == 0){
	print "File is empty !\n";exit;
}else{
  print "Analyse of the gff file $gff\n";
	foreach my $level (keys %{$hash_omniscient}){
    my $sizeLevel = keys %{$hash_omniscient->{$level}};
    print "$level we found $sizeLevel feature type:\n";
    foreach my $type (keys %{$hash_omniscient->{$level}}){
      print "$type\n";
      
      #useful to loop over level2 name
      if($level eq "level3"){
        foreach my $featureName (keys %{$hash_omniscient->{$level}{$type}}){
          $Level2featuresNames{$featureName}++;
        }

      }
    }
	}
}

### Handle to not print to much warning
my %WARNS;
my $nbWarn=5;
local $SIG{__WARN__} = sub {
  my $message = shift;
  my @thematic=split /@/,$message ;
   
  $WARNS{$thematic[0]}++;
  if ($WARNS{$thematic[0]} < $nbWarn){
    print $message;
  }
  elsif($WARNS{$thematic[0]} == $nbWarn){
    print "$thematic[0] ************** Too much WARNING message we skip the next **************\n";
  }
};


#######
# WE START FROM LEVEL3
# Should recreate level2 and level1 for exon-UTR-CDS
# Should create exon for CDS if not present (should take in account UTR)
######
print "CHECK STARTING FROM LEVEL3\n";
my %cptType;
my $newGeneIDcpt;
my @l3_endprocess;

foreach my $type (keys %{$hash_omniscient->{'level3'}}){
  foreach my $idL2 (keys %{$hash_omniscient->{'level3'}{$type}}){
    my $firstRound="yes";

    foreach my $f_l3 (sort {$a->start <=> $b->start}  @{$hash_omniscient->{'level3'}{$type}{$idL2}}) {
      
      my $mRNA_ID;
      my $ID_l3=$f_l3->_tag_value('ID');
      #################
      # Manage level3 #
      #################

      #################
      # Generate parent If needed
      if($f_l3->has_tag('Parent')){
        $mRNA_ID=$f_l3->_tag_value('Parent');
      }
      else{
        $cptType{$type}{$ID_l3}++;
        my $mrna_id_tmp=$ID_l3;
        $mrna_id_tmp=~s/-$type[-\.0-9]*//i;
        $mRNA_ID=$mrna_id_tmp."-mRNA".$cptType{$type}{$ID_l3};
        $f_l3->add_tag_value('Parent',$mRNA_ID);
        push(@{$hash_omniscient->{'level3'}{lc($type)}{lc($mRNA_ID)}}, $f_l3); #save feature with correct Parent
      }
      
      #### 
      # If first level3 feature we check if it's a CDS or exon and the counter part is not existing (exon or CDS) we create it.
      # /!\  $mRNA_ID and $idL2 are different, in case where ID doesn't exist in gff file (prokka case) we create  a undefined one that will not reuse then because we then give to it the correct ID decided. 
      my $L3_exists=undef;
      my $exon_exists=undef;
      my $cds_exists=undef;

      if ($firstRound){

        foreach my $type (keys %{$hash_omniscient->{'level3'}}){
          
        if(type_exist_in_omniscient($hash_omniscient, $type, 'level3',$mRNA_ID)){
            $L3_exists="yes";
          }
        }
        if(type_exist_in_omniscient($hash_omniscient, 'exon', 'level3',$mRNA_ID)){
            $exon_exists="yes";
        }
        if(type_exist_in_omniscient($hash_omniscient, 'cds', 'level3',$mRNA_ID)){
            $cds_exists="yes";
        }

        if($cds_exists and ! $exon_exists){$cptgg++;
          warn "WARNING create exon because we have CDS but no exon @ $cptgg\n";
          foreach my $cds_feature (@{$hash_omniscient->{'level3'}{'cds'}{$idL2}}){

            my $exon_feature = create_exon_from_cds($cds_feature,$mRNA_ID);
            push(@{$hash_omniscient->{'level3'}{'exon'}{lc($mRNA_ID)}}, $exon_feature);#save that exon
          }
        }
        $firstRound=undef;
      }
      
      #######################
      # Now create level 2 feature if needed
      if(! level2_exist_in_omniscient($hash_omniscient,$idL2)){
          
          # We save it to handle at the end and try to associate a parent
          #print "test L3_exists $L3_exists  exon_exists $exon_exists cds_exists $cds_exists\n";
          if(! $cds_exists and ! $exon_exists) {
            push(@l3_endprocess,$f_l3);
          }   
          else{ #We create parents

          my $mRNA_feature=clone($f_l3);#create a copy of the feature
          $mRNA_feature->primary_tag('mRNA');
          create_or_replace_tag($mRNA_feature,'ID', $mRNA_ID);

          # By default keep only ID and Parent
          if(! $keepAllAtt){
            remove_extra_attributes($mRNA_feature); 
          }
          # Handle attribute if asked (keep some, change tag name of them)
          if($attributes){
            manage_attributes($f_l3, $mRNA_feature);
          }


          ## link to level1 or create one if needed
          if(! overlap_existing_gene(\%hash_sortBySeq, \%$hash_omniscient, $mRNA_feature)){
              my $gene_feature=clone($mRNA_feature);#create a copy of the feature
              $gene_feature->remove_tag('Parent');
              $gene_feature->primary_tag('gene');
              my $gene_id="gene-".$newGeneIDcpt;
              create_or_replace_tag($gene_feature,'ID', $gene_id);
              $hash_omniscient->{'level1'}{'gene'}{$gene_id} = $gene_feature;#save that gene         

              #Link gene to the mRNA
              create_or_replace_tag($mRNA_feature,'Parent', $gene_id); #update mRNA
              push(@{$hash_omniscient->{'level2'}{'mrna'}{lc($gene_id)}}, $mRNA_feature);#save that mRNA

              $newGeneIDcpt++;
          }
        }
      } 
    }
  }
}

print "\nCHECK STARTING FROM LEVEL2\n";
#######
# WE START FROM LEVEL2
######
foreach my $type (keys %{$hash_omniscient->{'level2'}}){
  foreach my $idL1 (keys %{$hash_omniscient->{'level2'}{$type}}){
   
    my $firstRound="yes";

    foreach my $f_l2 (sort {$a->start <=> $b->start}  @{$hash_omniscient->{'level2'}{$type}{$idL1}}) {
      
      my $gene_ID;
      my $ID_l2=$f_l2->_tag_value('ID');
      
      #################
      # Generate parent
      if($f_l2->has_tag('Parent')){
        $gene_ID=$f_l2->_tag_value('Parent');
      }
      else{
        $cptType{$type}{$ID_l2}++;
        my $gene_id_tmp=$ID_l2;
        $gene_id_tmp=~s/-mrna[-\.0-9]*//i;
        $gene_ID=$gene_id_tmp."-gene".$cptType{$type}{$ID_l2};
        $f_l2->add_tag_value('Parent',$gene_ID);
      }
      
      #Manage level1
      if(lc($gene_ID) ne $idL1){  # In case where file didn't contain parent the both must be non equal
        if ($firstRound){         
            if(! type_exist_in_omniscient($hash_omniscient,'gene','level2', $gene_ID)){
             
              warn "WARNING create exon because we have $type without any\n";
              my $gene_feature=clone($f_l2);#create a copy of the feature
              $gene_feature->primary_tag('gene');
              $gene_feature->remove_tag('Parent');
              create_or_replace_tag($gene_feature,'ID',$gene_ID);

              # By default keep only ID and Parent
              if(! $keepAllAtt){
                remove_extra_attributes($gene_feature); 
              }
              # Handle attribute if asked (keep some, change tag name of them)
                if($attributes){
                manage_attributes($f_l2, $gene_feature);
              }

              $hash_omniscient->{'level1'}{'gene'}{lc($gene_ID)} = $gene_feature;#save that gene
            }
            $firstRound=undef;
          }
      #save the mRNA in the correct place      
      push(@{$hash_omniscient->{'level2'}{'mrna'}{lc($gene_ID)}}, $f_l2);#save that mRNA
      }

      #Manage level3
      my $L3_exists=undef;
      my $exon_exists=undef;
      my $cds_exists=undef;
      foreach my $type (keys %{$hash_omniscient->{'level3'}}){
        
      if(type_exist_in_omniscient($hash_omniscient, $type, 'level3',$ID_l2)){
          $L3_exists="yes";
        }
      }
      if(type_exist_in_omniscient($hash_omniscient, 'exon', 'level3',$ID_l2)){
          $exon_exists="yes";
      }
      if(type_exist_in_omniscient($hash_omniscient, 'cds', 'level3',$ID_l2)){
          $cds_exists="yes";
      }

      # CDS BUT NO EXON we create exons
      if($cds_exists and ! $exon_exists){
        foreach my $cds_feature (@{$hash_omniscient->{'level3'}{'cds'}{$ID_l2}}){

          my $exon_feature = create_exon_from_cds($cds_feature,$ID_l2);

          push(@{$hash_omniscient->{'level3'}{'exon'}{lc($ID_l2)}}, $exon_feature);#save that exon
        }
      }
      
      # No feature level3 exists !
      # We sould at least create one long exon corresponding to the mRNA
      if(! $L3_exists){ 

        my $exon_feature=clone($f_l2);#create a copy of the feature

        $exon_feature->primary_tag('exon');
        create_or_replace_tag($exon_feature,'Parent',$ID_l2);
        my $exon_id_tmp=$ID_l2;
        
        $exon_id_tmp=~s/-mrna[-\.0-9]*//i;
        $exon_id_tmp=$exon_id_tmp."-exon";

        create_or_replace_tag($exon_feature,'ID',$exon_id_tmp);

        # By default keep only ID and Parent
        if(! $keepAllAtt){
          remove_extra_attributes($exon_feature); 
        }
        # Handle attribute if asked (keep some, change tag name of them)
        if($attributes){
          manage_attributes($f_l2, $exon_feature);
        }

        push(@{$hash_omniscient->{'level3'}{'exon'}{lc($ID_l2)}}, $exon_feature);#save that exon

      }
    }
  }
}


### NOW check all L3 not linked to any feature if they can be

# sort by seq id
my %hash_sortBySeq_l2;
foreach my $tag_level2 (keys %{$hash_omniscient->{'level2'}}){
  foreach my $level2_id (keys %{$hash_omniscient->{'level2'}{$tag_level2}}){
    foreach my $featureL2 (@{$hash_omniscient->{'level2'}{$tag_level2}{$level2_id}}){
      my $position=$featureL2->seq_id."".$featureL2->strand;
      push (@{$hash_sortBySeq_l2{$tag_level2}{$position}}, $featureL2);
    }
  }
}

#add position
my %cptType;
foreach my $l3_feature (@l3_endprocess){
  my $overlap_something=undef;
  my $position_l3_feature=$l3_feature->seq_id."".$l3_feature->strand;

  foreach my $tag_l2 (keys %hash_sortBySeq_l2){
      
    if(exists ($hash_sortBySeq_l2{$tag_l2}{$position_l3_feature})){
      foreach my $featurel2 ( @{$hash_sortBySeq_l2{$tag_l2}{$position_l3_feature}}){

        if(($featurel2->start <= $l3_feature->end) and ($featurel2->end >= $l3_feature->start )){ # they overlap
          $overlap_something="yes";

          my $ID_l2=$featurel2->_tag_value('ID');
          
          my $new_l3_feature=clone($l3_feature);#create a copy of the feature

          create_or_replace_tag($new_l3_feature,'Parent',$ID_l2);
          my $l3_primay_tag=$new_l3_feature->primary_tag;
          my $ID_l3=$new_l3_feature->_tag_value('ID');

          $cptType{$l3_primay_tag}{$ID_l3}++;

          my $id_tmp=$ID_l3;
          $id_tmp=~s/-$l3_primay_tag[-\.0-9]*//i;
          $id_tmp=$id_tmp."-$l3_primay_tag".$cptType{$l3_primay_tag}{$ID_l3};
          
          create_or_replace_tag($new_l3_feature,'ID',$id_tmp);

          push(@{$hash_omniscient->{'level3'}{$l3_primay_tag}{lc($ID_l2)}}, $new_l3_feature);#save that exon
        }
      }
    }
  }
  if(! $overlap_something){
    warn "WARNING !! These feature overlap none of l2 feature ! @ ".$l3_feature->gff_string."\n"; exit;
  }
}
              


 print "END\n";

print_omniscient($hash_omniscient, $gffout);
$gffout->close();

print "well done - END\n";

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

sub create_exon_from_cds{

    my ($cds_feature, $mRNA_ID)=@_;

                my $exon_feature=clone($cds_feature);#create a copy of the feature

                $exon_feature->primary_tag('exon');
                create_or_replace_tag($exon_feature,'Parent',$mRNA_ID);
                my $exon_id_tmp=$exon_feature->_tag_value('ID');
                #print "exon_id_tmp1 $exon_id_tmp\n";
                $exon_id_tmp=~s/-cds/-exon/i;
                create_or_replace_tag($exon_feature,'ID',$exon_id_tmp);
                
                # By default keep only ID and Parent
                if(! $keepAllAtt){
                  remove_extra_attributes($exon_feature); 
                }
                # Handle attribute if asked (keep some, change tag name of them)
                if($attributes){
                  manage_attributes($cds_feature, $exon_feature);
                }

    return $exon_feature;
}


sub overlap_existing_gene{

  my($hash_sortBySeq, $hash_omniscient, $mRNA_feature)=@_;

  my $mrna_id_to_check = lc($mRNA_feature->_tag_value('ID')); 
  my ($start2,$end2) = get_longest_cds_start_end_from_mRNA_ID($hash_omniscient,$mrna_id_to_check); # look at CDS becaus ewe want only ioverlapinng CDS

  foreach my $seqid (keys %{$hash_sortBySeq}){
    
    foreach my $gene_feature ( @{$hash_sortBySeq{$seqid}}){
      
      my $gene_id = lc($gene_feature->_tag_value('ID')); 
      my ($start1,$end1) = get_longest_cds_start_end_from_geneID($hash_omniscient,$gene_id); # look at CDS because we want only ioverlapinng CDS

      
      if( ($start2 <= $end1) and ($end2 >= $start1) ){ #feature overlap considering extrem start and extrem stop. It's just to optimise the next step. Avoid to do the next step every time. So at the end, that test (current one) could be removed
           
        #now check at each CDS feature independently
        if (gene_mrna_features_overlap_by_cds($hash_omniscient,$gene_id, $mrna_id_to_check)){
          print "These two features overlap! :\n".$gene_feature->gff_string."\n".$mrna_id_to_check."\n";
          create_or_replace_tag($mRNA_feature,'Parent', lc($gene_id)); #update mrNA
          return 1;
        }
        else{return 0; }
      }
    }
  }
}

sub gene_mrna_features_overlap_by_cds{
  my  ($hash_omniscient,$gene_id, $mrna_id2)=@_;
  my $resu=undef;

  #check full CDS for each mRNA
  foreach my $mrna_feature (@{$hash_omniscient->{'level2'}{'mrna'}{$gene_id}}){

    my $mrna_id1 = lc($mrna_feature->_tag_value('ID'));     
   
    #check all cds pieces
    foreach my $cds_feature1 (@{$hash_omniscient->{'level3'}{'cds'}{$mrna_id1}}){
      foreach my $cds_feature2 (@{$hash_omniscient->{'level3'}{'cds'}{$mrna_id2}}){
          
        if(($cds_feature2->start <= $cds_feature1->end) and ($cds_feature2->end >= $cds_feature1->start )){ # they overlap
          $resu="yes";last;
        }
      }
      if($resu){last;}
    }
    if($resu){last;}
  }

  return $resu;
}

sub get_longest_cds_start_end_from_geneID{
  my  ($hash_omniscient,$gene_id)=@_;
  my $resu_start=100000000000;
  my $resu_end=0;

  #check full CDS for each mRNA
  foreach my $mrna_feature (@{$hash_omniscient->{'level2'}{'mrna'}{$gene_id}}){
    my $mrna_id = lc($mrna_feature->_tag_value('ID'));
    my $extrem_start=100000000000;
    my $extrem_end=0;
   
    #check all cds pieces
    foreach my $cds_feature (@{$hash_omniscient->{'level3'}{'cds'}{$mrna_id}}){
      if ($cds_feature->start < $extrem_start){
        $extrem_start=$cds_feature->start;
      }
      if($cds_feature->end > $extrem_end){
              $extrem_end=$cds_feature->end ;
      }
    }
    
    if($extrem_start < $resu_start){
        $resu_start=$extrem_start;
    }
    if($extrem_end > $resu_end){
      $resu_end=$extrem_end;
    }
  }
  return $resu_start,$resu_end;
}

sub get_longest_cds_start_end_from_mRNA_ID{
  my  ($hash_omniscient,$mrna_id)=@_;
  my $resu_start=100000000000;
  my $resu_end=0;

  #check full CDS for each mRNA
  my $extrem_start=100000000000;
   my $extrem_end=0;
   
  #check all cds pieces
  foreach my $cds_feature (@{$hash_omniscient->{'level3'}{'cds'}{$mrna_id}}){
      if ($cds_feature->start < $extrem_start){
        $extrem_start=$cds_feature->start;
      }
      if($cds_feature->end > $extrem_end){
        $extrem_end=$cds_feature->end ;
      }
  }

  return $extrem_start,$extrem_end;
}

sub level2_exist_in_omniscient{
  my ($hash_omniscient,$id_to_check)=@_;

  foreach my $type_l2 (keys %{$hash_omniscient->{'level2'}}){
    foreach my $id_l1 (keys %{$hash_omniscient->{'level2'}{$type_l2}}){
      if(exists($hash_omniscient->{'level2'}{$type_l2}{$id_l1})){
        foreach my $f_l2 (@{$hash_omniscient->{'level2'}{$type_l2}{$id_l1}}){
          if(lc($f_l2->_tag_value('ID')) eq lc($id_to_check)){
            return 1;
          }
        }
      }
    }   
  }
  return 0;
}

sub type_exist_in_omniscient {
  my ($hash_omniscient, $type, $level, $idL2)=@_;

  if(exists ($hash_omniscient->{$level})){
    if(exists ($hash_omniscient->{$level}{$type})){
      if(exists ($hash_omniscient->{$level}{$type}{lc($idL2)})){
        return 1;
      }
    }
  }

  return 0;
}


sub manage_attributes{

  my ($feature, $newfeat)=@_;

    #Now replace those asked
    foreach my $att (keys %attListOk){
      if ($feature->has_tag($att)){ 
         
        my @values=$feature->get_tag_values($att);
        my $value = shift @values ;
          
        if ($attListOk{$att} eq "null" ){ # the attribute name is kept inctact
          create_or_replace_tag($newfeat,$att, $value);
        }
        else{ # We replace the attribute name
          my $newAttributeName=$attListOk{$att};
          create_or_replace_tag($newfeat,$newAttributeName, $value);
        }
      } 
    }
}

sub remove_extra_attributes{
  my ($feature)=@_;

  my @list_tags=$feature->get_all_tags();

  foreach my $tag (@list_tags){
    if(lc($tag) ne "id" and lc($tag) ne "parent"){
      $feature->remove_tag($tag);
    }
  }
}

__END__

=head1 NAME

gff_2_full_gff.pl -
The script take a gff file as input, and will complete it if needed to create a full gff3 file.
If level3 or level2 or level1 feature are missing, we recontstruct them.
CDS whitout exon => create exon
Level2 wihtout exon => create exon
Level2 wihtout level1 => create level1
level3 whitout level2 => create level2
level3 whitout level1 => create level1
TESTED only on PROKKA gff that contains only CDS feature.
/!\ Deal with case where no PARENT attribute !

=head1 SYNOPSIS

    ./gff_2_full_gff.pl --gff=infile.gff [ -o outfile ]
    ./gff_2_full_gff.pl --help

=head1 OPTIONS

=over 8

=item B<-gff>, B<--infile>, B<--in>

Input gff file that will be processed.

=item B<--attributes>, B<--attribute>, B<--att>, B<-a>

Option usefull only if features from level1 or level2 have to be created to perform a proper gff3 file (By default only ID and PARENT are kept to create a correct gff3 file). 
In that case attributes specified, will be lift from the level down (level2 or level3) to the level up (level1 or level2).
/!\\ You must use "" if name contains spaces.
To replace the attribute name by a new attribute name you must use this formuation attributeName/newAttributeName.

=item B<--kaa> 
"Keep All Attribute", as the (abbreviation) name means, this option allow keeping all attribute. In that case, the --attribute(s) option is usefull ONLY in order to replace their tag names

=item B<-o> , B<--output> , B<--out>, B<--gff> or B<--outfile>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.

=back

=cut

