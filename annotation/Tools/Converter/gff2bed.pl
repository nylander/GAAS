#!/usr/bin/env perl

## BILS 2015 (www.bils.se)
## jacques.dainat@bils.se

use strict;
use Getopt::Long;
use BILS::Handler::GXFhandler qw(:Ok);
use BILS::Handler::GFF3handler qw(:Ok);
use Bio::Tools::GFF;
use Pod::Usage;

my $header = qq{
########################################################
# BILS 2015 - Sweden                                   #  
# jacques.dainat\@bils.se                               #
# Please cite BILS (www.bils.se) when using this tool. #
########################################################
};

my $outfile = undef;
my $gff = undef;
my $attributes = undef ;
my $help;

if( !GetOptions(
    "help" => \$help,
    "gff=s" => \$gff,
    "outfile=s" => \$outfile))
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
           -message => "$header\nAt least 1 parameters is mandatory.\n",
           -verbose => 0,
           -exitval => 1 } );
}

## Manage output file
my $bedout;
if ($outfile) {
open(my $fh, '>', $outfile) or die "Could not open file '$outfile' $!";
  #$bedout= Bio::FeatureIO->new(-fh => $fh, -format => 'bed' );
  $bedout=$fh;
}
else{
  #$bedout = Bio::FeatureIO->new(-fh => \*STDOUT,  -format => 'bed');
  $bedout=\*STDOUT ;
}

### Parse GTF input file 
my ($hash_omniscient, $hash_mRNAGeneLink) = BILS::Handler::GFF3handler->slurp_gff3_file_JD($gff);
# END parsing


# NOT USED BECAUSE in Bio::FeatureIO::bed:
#  my $block_count = '';  #not implemented, used for sub features
#  my $block_sizes = '';  #not implemented, used for sub features
#  my $block_starts = ''; #not implemented, used for sub feature
#   #################
#   # == LEVEL 1 == #
#   #################
#   foreach my $tag_l1 (keys %{$hash_omniscient->{'level1'}}){ # tag_l1 = gene or repeat etc...
#     foreach my $id_l1 (keys %{$hash_omniscient->{'level1'}{$tag_l1}}){

#       #my feature is a Bio::SeqFeature::Generic
#       my $feature_l1=$hash_omniscient->{'level1'}{$tag_l1}{$id_l1};

#       #create a new  Bio::SeqFeature::Annotated object;
#       my $newObj = Bio::SeqFeature::Annotated->new();
#       #initialize this object with the contents of another feature
#       $newObj->from_feature($feature_l1);
#       #print the new object
#       print $bedout->write_feature($newObj);
#   }
# }


#################
# == LEVEL 1 == #
#################
foreach my $tag_l1 (keys %{$hash_omniscient->{'level1'}}){ # tag_l1 = gene or repeat etc...
  foreach my $id_l1 (keys %{$hash_omniscient->{'level1'}{$tag_l1}}){

    my $feature_l1=$hash_omniscient->{'level1'}{$tag_l1}{$id_l1};
    my $size_l1 = $feature_l1->end - $feature_l1->start;

    print $bedout $feature_l1->seq_id."\t".$feature_l1->start."\t".$feature_l1->end."\t".$feature_l1->primary_tag."\t".$feature_l1->score."\t".$feature_l1->strand."\t".$feature_l1->start."\t".$feature_l1->end."\t0\t1\t".$size_l1."\t0\n";

      #################
      # == LEVEL 2 == #
      #################
      foreach my $tag_l2 (keys %{$hash_omniscient->{'level2'}}){ # tag_l2 = mrna or mirna or ncrna or trna etc...         
        if ( exists ($hash_omniscient->{'level2'}{$tag_l2}{$id_l1} ) ){  
          foreach my $feature_l2 ( @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}}) {
            my $size_l2 = $feature_l2->end - $feature_l2->start;
            print $bedout $feature_l2->seq_id."\t".$feature_l2->start."\t".$feature_l2->end."\t".$feature_l2->primary_tag."\t".$feature_l2->score."\t".$feature_l2->strand."\t".$feature_l2->start."\t".$feature_l2->end."\t0\t1\t".$size_l2."\t0\n";

            #################
            # == LEVEL 3 == #
            #################
            my $level2_ID = lc( $feature_l2->_tag_value('ID') );

            foreach my $tag_l3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
              if ( exists ($hash_omniscient->{'level3'}{$tag_l3}{$level2_ID} ) ){

                my $first_feature_l3;
                my $last_feature_l3;
                my $final_sizeList; 
                my $final_startList;

                my $originStart;
                my $cpt = 0 ;
                my $nb_feat = 0 ;
                my $score = 0 ;

                foreach my $feature_l3 ( sort { $a->start <=> $b->start } @{$hash_omniscient->{'level3'}{$tag_l3}{$level2_ID}}) {
                  
                  my $size_l3 = $feature_l3->end - $feature_l3->start;
                  
                  
                  if($cpt==0){ 
                    $first_feature_l3 = $feature_l3;
                    $originStart = $feature_l3->start;                   

                    $final_sizeList.="$size_l3";
                    $final_startList.="0";
                    $cpt++;
                  }
                  else{
                      my $start_corrected= $feature_l3->start - $originStart;
                      $final_sizeList.=",$size_l3";
                      $final_startList.=",$start_corrected";
                  }
                  $last_feature_l3 = $feature_l3;
                  $nb_feat++; 
                  $score += $feature_l3->score;
                }
                $score = $score/$nb_feat;
                print $bedout $first_feature_l3->seq_id."\t".$first_feature_l3->start."\t".$last_feature_l3->end."\t".$first_feature_l3->primary_tag."\t".$score."\t".$first_feature_l3->strand."\t".$first_feature_l3->start."\t".$last_feature_l3->end."\t0\t1\t".$final_sizeList."\t".$final_startList."\n";
            }
          }
        }
      }
    }
  } 
}

__END__

