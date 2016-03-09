#!/usr/bin/perl

## TO DO:
## Dont manage if some isoform (not all) should be modelate for a same gene. (gene will appear 2 times)
# If tehre is only the level2 to build

use strict;
use Getopt::Long;
use Clone 'clone';
use Pod::Usage;
use Bio::Tools::GFF;
use BILS::Handler::GTFhandler qw(:Ok);
use BILS::Handler::GXFhandler qw(:Ok);
use BILS::Handler::GFF3handler qw(:Ok);
use Data::Dumper;

my $header = qq{
########################################################
# BILS 2015 - Sweden                                   #  
# jacques.dainat\@bils.se                               #
# Please cite BILS (www.bils.se) when using this tool. #
########################################################
};

my $outfile = undef;
my $gtf = undef;
my $attributes=undef;
my $keepAllAtt=undef;
my $help= 0;

if ( !GetOptions(
    "help|h" => \$help,
    "gtf|in|infile=s" => \$gtf,
    "attributes|a|att=s" => \$attributes,
    "kaa" => \$keepAllAtt,
    "outfile|out|o|output|gff=s" => \$outfile))

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
 
if ( ! (defined($gtf)) ){
    pod2usage( {
           -message => "\nAt least 1 parameter is mandatory:\nInput reference gtf file (--f)\n",
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
  @attListPair= split(/,/, $attributes);
  my $nbAtt=$#attListPair+1;
  print "In addition of gene_id attribute, we will keep $nbAtt attribute(s):\n";
  foreach my $attributeTuple (@attListPair){ 
    my @attList= split(/\//, $attributeTuple);
    if($#attList == 0){ # Attribute alone
      $attListOk{$attList[0]}="null";
      print "$attList[0]\n";
    }
    else{ # Attribute we have to replace by a new name
      $attListOk{$attList[0]}=$attList[1];
      print "$attList[0] (Will be replaced by $attList[1])\n";
    }
  }
  print "\n";
}



### Parse GTF input file 
my ($hash_omniscient, $hash_mRNAGeneLink) = BILS::Handler::GTFhandler->slurp_gtf_file_JD($gtf);

## Check duplication
#my $size = keys %$hash_duplicatedFeature;
#if($size == 0){
#	print "No duplicated feature found !\n";
#}else { print "//// !! ACHTUNG \\\\\\\\\ $size duplicated feature found !";}

my $hash_level1=$hash_omniscient->{"level1"};
my $hash_level2=$hash_omniscient->{"level2"};
my $hash_level3=$hash_omniscient->{"level3"};
my %Level2featuresNames;

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
  print "Analyse of the gtf file $gtf\n";
	foreach my $level (keys %$hash_omniscient){
    my @keysLevel = keys (%{$hash_omniscient->{$level}});
    my $sizeLevel = @keysLevel;
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

if($sizeLevel1 == 0 and $sizeLevel2 == 0){
  print "Level1 and Leve2 featured must be build for all features!!\n\n";
}
#my %hash=%$hash_omniscient;
#print Dumper %hash;
#__END__
### check omniscient


###########
# CONSTANTS
my $NewGeneID=0;

#######
# WE START FROM LEVEL2
######
my %level1_alreadyAnnotated;
my @meta_listfeature;
foreach my $parent (keys %Level2featuresNames){


  ################################
  # CHECK if level2 features exist
  if (exists($hash_mRNAGeneLink->{$parent})){ #parent attribute exists for level2 features
    
    my $gene_id=$hash_mRNAGeneLink->{$parent};
    my @gene_id_list=($gene_id);

    ############
    # CHECK if level1 features exist else we will build it
    if (exists($hash_omniscient->{'level1'})){ # do we have level1 features ?
      foreach my $tag (keys %{$hash_omniscient->{'level1'}}){ 
        if (exists($hash_omniscient->{'level1'}{$tag}{$gene_id})){ #check now if parent really exists as feature

          if(! exists ($level1_alreadyAnnotated{$hash_mRNAGeneLink->{$parent}}) ){
            
            ######
            # Now print it
            gtf2gff_features_in_omniscient_from_level1_id_list($hash_omniscient, \@gene_id_list, $gffout);
            print_omniscient_from_level1_id_list($hash_omniscient, \@gene_id_list, $gffout);
            $level1_alreadyAnnotated{$hash_mRNAGeneLink->{$parent}}++;
          }
        }
        ########
        # BUILD LEVEL1 feature
        # case where level1 exits but not the feature 
        else{
          foreach my $tag (keys %{$hash_omniscient->{'level2'}}){ 
            if (exists ($hash_omniscient->{'level2'}{$tag}{$gene_id})){
              my $level2_feature=@{$hash_omniscient->{'level2'}{$tag}{$gene_id}}[0];
              my $level1_feature=clone($level2_feature);#create a copy of one of the level2 feature
              my $tag_level1='gene';
              $level1_feature->primary_tag($tag_level1); # change promary tag. Assume is only gene feature !!!!
              $hash_omniscient->{"level1"}{$tag_level1}{$gene_id}=$level1_feature; # push the feature in omniscient now!

              #######
              # Manage attribute

              # By default keep only ID and Parent
              if($keepAllAtt){
                lift_all_attributes($level2_feature, $level1_feature); 
              }

              if($attributes){
                manage_attributesa($level2_feature, $level1_feature);
              }

              ######
              # Now print it
              gtf2gff_features_in_omniscient_from_level1_id_list($hash_omniscient, \@gene_id_list, $gffout);
              print_omniscient_from_level1_id_list($hash_omniscient, \@gene_id_list, $gffout);
              $level1_alreadyAnnotated{$hash_mRNAGeneLink->{$parent}}++;
            }      
          }
        }
      }
    }
    ########
    # BUILD LEVEL1 feature
    # case where level1 still not exits !!
    else{ # we have to create the gtf feature for gene
      foreach my $tag (keys %{$hash_omniscient->{'level2'}}){ 
        if (exists ($hash_omniscient->{'level2'}{$tag}{$gene_id})){
          my $level2_feature = @{$hash_omniscient->{'level2'}{$tag}{$gene_id}}[0];
          my $level1_feature=clone($level2_feature);#create a copy of one of the level2 feature
          my $tag_level1='gene';
          $level1_feature->primary_tag($tag_level1); # change promary tag. Assume is only gene feature !!!!
          $hash_omniscient->{"level1"}{$tag_level1}{$gene_id}=$level1_feature; # push the feature in omniscient now!

          #######
          # Manage attribute

          # By default keep only ID and Parent
          if($keepAllAtt){
            lift_all_attributes($level2_feature, $level1_feature); 
          }

          if($attributes){
             manage_attributesa($level2_feature, $level1_feature);
          }

          ######
          # Now print it
          gtf2gff_features_in_omniscient_from_level1_id_list($hash_omniscient, \@gene_id_list, $gffout);
          print_omniscient_from_level1_id_list($hash_omniscient, \@gene_id_list, $gffout);
          $level1_alreadyAnnotated{$hash_mRNAGeneLink->{$parent}}++;
        }      
      }
    }
  }
  

  ##########
  # WE have to reconstruc level2 and level1
  #########
  else{
    my @L3_f_List=();
    my $feature_withFunc;

    #get level3 feature list
    foreach my $tag (keys %{$hash_level3}){
#     
        if(exists ($hash_level3->{$tag}{$parent}) ){
          @L3_f_List = (@L3_f_List, @{$hash_level3->{$tag}{$parent}});

          # save one feature to lift Attributes
          if($tag eq "exon" or $tag eq "cds"){
            $feature_withFunc = @{$hash_level3->{$tag}{$parent}}[0];
          }
        }
    }
    if(! @L3_f_List)
    { print "No level3 feature present for this one: $parent\n";next;}

    my ($ref_list_gene, $ref_list_mRNA, $NewGeneIDtmp) = reconstruct_locus_without_transcripts_with_seq_id(\@L3_f_List, $NewGeneID);
    $NewGeneID=$NewGeneIDtmp; # keep track If Gene Id modified

    ######
    # Manage attribute

    # By default keep only ID and Parent
    if($keepAllAtt){    
      foreach my $feature (@$ref_list_gene){
        lift_all_attributes($feature_withFunc,$feature);
      }
      foreach my $feature (@$ref_list_mRNA){
        lift_all_attributes($feature_withFunc,$feature);
      }
    }

    if($attributes){
      foreach my $feature (@$ref_list_gene){
        manage_attributes($feature_withFunc, $feature);
      }
      foreach my $feature (@$ref_list_mRNA){
        manage_attributes($feature_withFunc, $feature);
      }
    }

    my @ref_list_Level1Level2_features;
    # check if we have buit gene feature several times
    my $gene_id = $ref_list_gene->[0]->_tag_value('ID');

    if(! exists ($level1_alreadyAnnotated{$gene_id}) ){
      $level1_alreadyAnnotated{$gene_id}++;
      @ref_list_Level1Level2_features=(@$ref_list_gene, @$ref_list_mRNA);
    }
    else{
      # we dont use the gene feature because it already exist
      @ref_list_Level1Level2_features=(@$ref_list_mRNA);
    }

    # Now we have level1 and level2 ready

    ###########
    # Manage level3 features
    ###########
    my @list_level3_feature;
    foreach my $type (keys %{$hash_level3}){
      if( exists($hash_level3->{$type}{$parent}) ){
        push(@list_level3_feature, @{$hash_level3->{$type}{$parent}});   
      }
    }
    my $newlist_level3_feature = manage_level3_attributes_ID_parent(@list_level3_feature);

    #PRINT all the bucket
    push (@meta_listfeature, (@ref_list_Level1Level2_features, @$newlist_level3_feature));
	}

}

#If some has been reconstructed - print them
if ($#meta_listfeature != -1){ 
  my ($hash_omniscient_sorted, $hash_mRNAGeneLink_sorted) = create_omniscient_from_feature_list(\@meta_listfeature);
  print_omniscient($hash_omniscient_sorted, $gffout);
}

$gffout->close();

print "well done - END\n";


sub manage_attributes{

  my ($feature, $newfeat)=@_;

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

sub lift_all_attributes{
  my ($feature, $newfeat)=@_;

  my @list_tags=$feature->get_all_tags();

  foreach my $tag (@list_tags){
    if(lc($tag) ne "id" and lc($tag) ne "parent"){
      my @values=$feature->get_tag_values($tag);
      create_or_replace_tag($newfeat,$tag, @values);
    }
  }
}

__END__

=head1 NAME

gtf2gff_universal.pl -
The script take a gtf file as input, and will translate it in gff3.
If level2 feature are missing, we assume that level1 feature are also missing, so we recontruct them (old ensembl gtf format).

=head1 SYNOPSIS

    ./gtf2gff3_universal.pl -gtf=infile.gff [ -o outfile ]
    ./gtf2gff3_universal.pl --help

=head1 OPTIONS

=over 8

=item B<-gtf>, B<--infile>, B<--in>

Input gtf file that will be convert.

=item B<-attributes>, B<--att>, B<-a>

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

