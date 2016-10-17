#!/usr/bin/env perl

## BILS 2015
## jacques.dainat@bils.se

use strict;
use Pod::Usage;
use Getopt::Long;
use Bio::SeqFeature::Generic;
use Bio::Tools::GFF;

my $usage = qq{
########################################################
# BILS 2015 - Sweden                                   #  
# jacques.dainat\@bils.se                               #
# Please cite BILS (www.bils.se) when using this tool. #
########################################################

Usage: perl my_script.pl --gtf Infile [--out outfile]
  Getting help:
    [--help]

  Input:
    [--gtf filename]
    The name of the gtf file to convert. 
    
  Ouput:    
    [--out filename]
        The name of the output file (A GFF file).

  Optional
  [--clean]
        Clean features a maximun

  [--attributes "attribute 1",attribute2,attribute3]
    The name of the attributes from input file that must be kept in the output file. By default only ID and PARENT are kept to create a correct gff3 file.
    /!\\ You must use "" if name contains spaces.
    To replace the attribute name by a new attribute name you must use this formuation attributeName/newAttributeName.


  /!\\  We refer to gtf3 because it's a gtf that contains all features. i.e No need to reconstruct transcript or gene from exons /!\\
};

my @EnsemblAttributeArray=("gene_id","transcript_id","exon_id","exon_number","exon_version","gene_biotype","gene_name","gene_source","gene_version","transcript_biotype","transcript_name","transcript_source","transcript_version");
my $outfile = undef;
my $gtf = undef;
my $clean = undef;
my $attributes = undef ;
my $help;

if( !GetOptions(
    "help" => \$help,
    "gtf=s" => \$gtf,
    "clean|c" => \$clean,
    "attributes|a|att=s" => \$attributes,
    "outfile=s" => \$outfile))
{
    pod2usage( { -message => "Failed to parse command line\n$usage",
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if ( ! (defined($gtf)) ){
    pod2usage( {
           -message => "Missing the --gtf argument\n$usage",
           -verbose => 0,
           -exitval => 2 } );
}

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

## Manage output file
my $gffout;
if ($outfile) {
open(my $fh, '>', $outfile) or die "Could not open file '$outfile' $!";
  $gffout= Bio::Tools::GFF->new(-fh => $fh, -gff_version => 3 );
}
else{
  $gffout = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 3);
}

### Read gtf input file /!\GTF, is also known as GFF v2.5. 
my $gtfio = Bio::Tools::GFF->new(-file => $gtf, -gff_version => 2.5);

# write new format
my $featureID=1;
while( my $feature = $gtfio->next_feature()) {
  my $newfeat = Bio::SeqFeature::Generic->new(
            -seq_id       => $feature->seq_id(),
            -start        => $feature->start(),
            -end          => $feature->end() ,
            -strand       => $feature->strand() , 
            -primary      => $feature->primary_tag() , # -primary_tag is a synonym
            -source_tag   => $feature->source_tag() ,
            -display_name => $feature->display_name() ,
            -score        => $feature->score() );

  # add attributes selected from the old features in new features
  if ($attributes){
    if (defined($clean)){ # IF THAT OPTION IS SET UP WE WILL REMOVE NON USEFUL INFORMATION (I.E GENE OR TRANSCRIPT INFORMATION IN EXON OR CDS FEATURES)
      foreach my $att (keys %attListOk){
        if( !((($att eq "transcript_id") and ($feature->primary_tag() eq "gene")))){ # particular case.  We never want to keep transcript information in the gene feature. Indeed a gene feature is allow to have several transcript.
          # we dont want to keep the gene_name transcript_name transcript_biotype and
        if( !((($att eq "gene_name") or ($att eq "gene_id") or ($att eq "gene_biotype")) and (($feature->primary_tag() eq "transcript") or ($feature->primary_tag() eq "mrna")))){ 
        if( !((grep {$_ eq $att} @EnsemblAttributeArray) and !((lc($feature->primary_tag()) eq "transcript") or (lc($feature->primary_tag()) eq "mrna") or (lc($feature->primary_tag()) eq "gene")))){ # we remove all information of features different of trranscript and gene
          if ($feature->has_tag($att)){
            my @values=$feature->get_tag_values($att);
            my $value = shift @values ;
            if ($attListOk{$att} eq "null" ){ # the attribute name is kept inctact
              $newfeat->add_tag_value($att,$value);
            }
            else{ # We replace the attribute name
              my $newAttributeName=$attListOk{$att};
              $newfeat->add_tag_value($newAttributeName,$value);
            }
          }
        }
        }  
        }
      }
    }
    else{
      foreach my $att (keys %attListOk){
        if ($feature->has_tag($att)){
          my @values=$feature->get_tag_values($att);
          my $value = shift @values ;
          if ($attListOk{$att} eq "null" ){ # the attribute name is kept inctact
            $newfeat->add_tag_value($att,$value);
          }
          else{ # We replace the attribute name
            my $newAttributeName=$attListOk{$att};
            $newfeat->add_tag_value($newAttributeName,$value);
          }
        } 
      }
    }
  }

  if ($feature->primary_tag() eq "gene"){
    my @ID=$feature->get_tag_values("gene_id");
    $newfeat->add_tag_value("ID", $ID[0]);
  }
  elsif($feature->primary_tag() eq "transcript"){
    my @ID=$feature->get_tag_values("transcript_id");
    my @parentID=$feature->get_tag_values("gene_id");
    $newfeat->add_tag_value("ID", $ID[0]);
    $newfeat->add_tag_value("Parent", $parentID[0]);
  }
  else{
    my @parentID=$feature->get_tag_values("transcript_id");
    $newfeat->add_tag_value("ID", $featureID);
    $newfeat->add_tag_value("Parent", $parentID[0]);
    $featureID++;
  }
  

	$gffout->write_feature($newfeat);		
}
$gtfio->close();
$gffout->close();
