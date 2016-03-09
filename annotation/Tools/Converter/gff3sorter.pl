#!/usr/bin/perl

## BILS 2015
## jacques.dainat@bils.se
use Try::Tiny;
use Pod::Usage;
use strict;
use Getopt::Long;
use Bio::SeqFeature::Generic;
use Bio::Tools::GFF;
use Pod::Usage;

my $usage = qq{
########################################################
# BILS 2015 - Sweden                                   #  
# jacques.dainat\@bils.se                               #
# Please cite BILS (www.bils.se) when using this tool. #
########################################################

Usage: perl my_script.pl --gff Infile [--out outfile]
  Getting help:
    [--help]

  Input:
    [--gff filename]
    The name of the gtf file to convert. 
    
  Ouput:    
    [--out filename]
        The name of the output file (A GFF file).

  Optional
  [--attributes "attribute 1",attribute2,attribute3]
    The name of the attributes from input file that must be kept in the output file. By default only ID and PARENT are kept to create a correct gff3 file.
    /!\\ You msut use "" if name contains spaces.
    To replace the attribute name by a new attribute name you must use this formulation attributeName/newAttributeName.


  /!\\  We refer to gtf3 because it's a gtf that contains all features. i.e No need to reconstruct transcript or gene from exons /!\\
};

my $outfile = undef;
my $gtf = undef;
my $attributes = undef ;
my $help;

if ( !GetOptions(
    "help" => \$help,
    "gff=s" => \$gtf,
    "outfile|out|o=s" => \$outfile))
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
           -message => "Missing the --gff argument\n$usage",
           -verbose => 0,
           -exitval => 2 } );
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
my $gtfio = Bio::Tools::GFF->new(-file => $gtf, -gff_version => 3);

my %hash_parent;
my %hash_children;
my %hash_link_parent_children;

# write new format
while( my $feature = $gtfio->next_feature()) {
    my @ID=$feature->get_tag_values("ID");
    my @tags = $feature->get_all_tags();
    my $ParentPresent="no";
    foreach my $tag (@tags){
      if($tag eq "Parent" ){
        $ParentPresent="yes";
      }
    }
    if($ParentPresent eq "yes"){
      my @Parent=$feature->get_tag_values("Parent");
      $hash_children{$ID[0]}=$feature;
      if(exists ( $hash_link_parent_children{$Parent[0]} ) ){
        push( $hash_link_parent_children{$Parent[0]},$ID[0] );

        my @size= @{$hash_link_parent_children{$Parent[0]}};
      }
      else{
        $hash_link_parent_children{$Parent[0]} =[$ID[0]];
        my @size= @{$hash_link_parent_children{$Parent[0]}};
      }
    }
    else{
      $hash_parent{$ID[0]}=$feature;
    }
}

# Print genes
foreach my $key (keys %hash_parent){
  $gffout->write_feature($hash_parent{$key});
  my @listchildren=@{$hash_link_parent_children{$key}};
  #Print transcripts
  foreach my $idChildren (@listchildren){ 
    $gffout->write_feature($hash_children{$idChildren});
    my @ID=$hash_children{$idChildren}->get_tag_values("ID");
    my @list_sub_children=@{$hash_link_parent_children{$ID[0]}};
    #Print other feature
    foreach my $id_sub_Children (@list_sub_children){
       $gffout->write_feature($hash_children{$id_sub_Children});
    }
  }
}
$gtfio->close();
$gffout->close();
