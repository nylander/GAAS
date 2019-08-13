#!/usr/bin/env perl


use Carp;
use strict;
use Clone 'clone';
use Getopt::Long;
use Data::Dumper;
use Pod::Usage;
use List::MoreUtils qw(uniq);
use Bio::Tools::GFF;
use BILS::Handler::GFF3handler qw(:Ok);
use BILS::Handler::GXFhandler qw(:Ok);

my $header = qq{
########################################################
# BILS 2018 - Sweden                                   #  
# jacques.dainat\@bils.se                               #
# Please cite BILS (www.bils.se) when using this tool. #
########################################################
};

my $outfile = undef;
my $gff = undef;
my $help= 0;

if ( !GetOptions(
    "help|h" => \$help,
    "gff|g=s" => \$gff,
    "output|outfile|out|o=s" => \$outfile))

{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($help) {
    pod2usage( { -verbose => 2,
                 -exitval => 2,
                 -message => "$header\n" } );
}

if (! $gff){
    pod2usage( {
           -message => "\nAt least 1 files is mandatory:\n --gff file1\n\n",
           -verbose => 0,
           -exitval => 2 } );
}

######################
# Manage output file #
my $gffout;
if ($outfile) {
  $outfile=~ s/.gff//g;
open(my $fh, '>', $outfile.".gff") or die "Could not open file '$outfile' $!";
  $gffout= Bio::Tools::GFF->new(-fh => $fh, -gff_version => 3 );
}
else{
  $gffout = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 3);
}


                                                      #######################
                                                      #        MAIN         #
#                     >>>>>>>>>>>>>>>>>>>>>>>>>       #######################       <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


my $verbose=1;
my $resume_case=undef;

######################
### Parse GFF input #
my ($omniscient, $mRNAGeneLink) = slurp_gff3_file_JD({ input => $gff
                                                              });  
print ("$gff GFF3 file parsed\n");
info_omniscient($omniscient);


# Quick stat hash before complement
my %quick_stat1;
foreach my $level ( ('level1', 'level2') ){
  foreach  my $tag (keys %{$omniscient->{$level}}) {
    my $nb_tag = keys %{$omniscient->{$level}{$tag}};
    $quick_stat1{$level}{$tag} = $nb_tag;
  }
}

########
# Sort the genes to loop over them from the left to right.
my $sortBySeq = gather_and_sort_l1_by_seq_id_and_strand($omniscient);
my %alreadyChecked;

foreach my $locusID ( keys %{$sortBySeq}){ # tag_l1 = gene or repeat etc...
  if ( exists_keys( $sortBySeq, ( $locusID ) ) ){
    foreach my $tag_l1 ( keys %{$sortBySeq->{$locusID}} ) { 
      my $sorted_locations = _create_sorted_list_of_locations(\@{$sortBySeq->{$locusID}{$tag_l1}});

      # Go through location from left to right ### !!
      for (my $i = 0; $i < scalar @$sorted_locations; $i++) {

        my $location = shift @$sorted_locations;# This location will be updated on the fly
        my $id_l1_left = $location->[0];
        if( l1_has_l3_type($omniscient, $id_l1_left, 'cds')){

          # Go through location from left to right ### !!
          foreach my $locationB ( @{$sorted_locations} ) {
            my $id_l1_right = $locationB->[0];
            if( l1_has_l3_type($omniscient, $id_l1_right, 'cds')){

              #If location_to_check start if over the end of the M1erence location, we stop
              if($locationB->[1] > $location->[2]) {last;} 

              # Let's check if neighbor gene overlap at Gene LEVEL !
              if (($location->[1] <= $locationB->[2]) and ($location->[2] >= $locationB->[1])){
                ### OVERLAP ###
                #let's check at UTR level
                _check_gene_overlap_at_UTR($omniscient , $id_l1_left, $id_l1_right, $verbose); #If contains UTR
              }
            }
          }
        }
      }
    }
  }
}
  print "We fixed $resume_case case where feature has been merged within the same locus\n" if($verbose >= 1 and $resume_case);

########
# Print results
print_omniscient($omniscient, $gffout);  

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

sub _create_sorted_list_of_locations{
  my  ($list)=@_;
  my @result;
  foreach my $feature (@{$list}){
    push(@result, [$feature->_tag_value('ID'), $feature->start, $feature->end ])
  }

  @result = sort { $a->[1] <=> $b->[1] }  @result;
  return \@result;
}


sub _check_gene_overlap_at_UTR{
  my  ($hash_omniscient, $gene_id, $gene_id2, $verbose)=@_;

  my $verbose=1;
  
  #collect extrem cds positions.
  my ($gene1_cds_start, $gene1_cds_end) = get_most_right_left_cds_positions($omniscient, $gene_id);
  my ($gene2_cds_start, $gene2_cds_end) = get_most_right_left_cds_positions($omniscient, $gene_id2);

  my $utr_mrna1 = undef;
  my $utr_mrna2 = undef;
  print "lets go: $gene_id, $gene_id2 \n";


  foreach my $l2_type (keys %{$hash_omniscient->{'level2'}} ){
    if(exists_keys($hash_omniscient,('level2', $l2_type, $gene_id))){
      foreach my $mrna_feature (@{$hash_omniscient->{'level2'}{$l2_type}{$gene_id}}){       
        my $righ_utrs = undef;
        my $sorted_cds = get_cds_from_l2($hash_omniscient, $mrna_feature);
        if ($sorted_cds){ #this l2 has a cds
          $righ_utrs = _get_right_utrs($hash_omniscient, $mrna_feature, $gene1_cds_end);
        }


        foreach my $l2_type_B (keys %{$hash_omniscient->{'level2'}} ){
          if(exists_keys($hash_omniscient,('level2', $l2_type_B, $gene_id2))){
            foreach my $mrna_feature_B (@{$hash_omniscient->{'level2'}{$l2_type}{$gene_id2}}){       
              my $sorted_cds_B = get_cds_from_l2($hash_omniscient, $mrna_feature_B);
              if ($sorted_cds_B){ #this l2 has a cds
                if($sorted_cds_B[0]->start <=  $sorted_cds[$#sorted_cds]->end and $sorted_cds_B[$#sorted_cds_B]->end  >= $sorted_cds[0]->start){
                  print "Both CDS overlap, we will not touch their UTRs\n";
                }
                my $left_utrs = _get_left_utrs($hash_omniscient, $mrna_feature_B, $gene2_cds_start);
                if($righ_utrs and $left_utrs){
                  #-----HERE BOTH HAVE UTR-----
                  my $length_utr_M1 =  $righ_utrs->[$#righ_utrs]->end - $righ_utrs->[0]->start +1;
                  my $length_utr_M2 =  $left_utrs->[$#left_utrs]->end  - $left_utrs->[0]->start + 1;
                  my $separting_point = undef;
                  if($length_utr_M1 > $length_utr_M2){
                    $separting_point = $M2_utr_left_start;
                  }
                  #print "($length_utr_M1 > $dist_between_M1_and_M2 and $length_utr_M2 > $dist_between_M1_and_M2)\n" if $verbose;
                  # If $dist_between_M1_and_M2 = 1 we cannot share one nucleotide between two utr. Both should have it
                  if($length_utr_M1 > $dist_between_M1_and_M2 and $length_utr_M2 > $dist_between_M1_and_M2 and $dist_between_M1_and_M2 > 1) {
                    $separting_point =  $M1_cds_end  + int($dist_between_M1_and_M2 / 2) + 1;
                    #print "separting_point  $separting_point $M2_cds_start - $M1_cds_end - 1 =  $dist_between_M1_and_M2 ".int($dist_between_M1_and_M2 / 2)."\n" if $verbose;
                  }
                  _shrink_utr_right($separting_point)
                  _shrink_utr_left($separting_point)
                }
                elsif($righ_utrs){
                  #-----HERE left gene has UTR to his right to check
                  my $separting_point = $sorted_cds_B->[0]->start;
                  _shrink_utr_right($separting_point)
                }
                elsif($left_utrs){
                  #-----HERE right gene has UTR to his left to check
                  my $separting_point = $sorted_cds->[$#$sorted_cds]->end;
                  _shrink_utr_left($separting_point)
                }
                else{
                  print "None of the two features have CDS" if ($verbose);
                }
              }
            }
          }
        }
      }
    }
  }
}




############################################################
#check if UTR right of model left is overlaping, and fix it.
sub _control_utr_from_model_left{
  my ($l2_feature_M1, $M1_cds_start, $M1_cds_end, $l2_feature_M2, $M2_cds_start, $M2_cds_end, $omniscient, $verbose)=@_;

  ####################################
  #peculiar case CDS1 and CD2 overlap. We dont touch that case
  if ($M1_cds_start <= $M2_cds_end and $M1_cds_end >= $M2_cds_start){
    print "peculiar case CDS1 and CD2 overlap. We dont touch that case\n";
    return;
  }

  my $l2_M1_id  = lc($l2_feature_M1->_tag_value('ID'));
  my $l2_M2_id  = lc($l2_feature_M2->_tag_value('ID')); 

  ##############################
  # Look at utr of the M1erence
  my ($M1_utr_left_start, $M1_utr_left_end, $M1_utr_right_start, $M1_utr_right_end) = _get_utrs_extremities($l2_M1_id, $M1_cds_end, $omniscient);
  if(! $M1_utr_right_start){return;} #No UTR at right of model1, nothing to modify

  ###########################
  # Look at utr of the M2
  my ($M2_utr_left_start, $M2_utr_left_end, $M2_utr_right_start, $M2_utr_right_end) = _get_utrs_extremities($l2_M2_id, $M2_cds_end, $omniscient);


  ##############################
  # CHECK RIGHT UTR OF MODEL 1 #
  ##############################
  my $utr_redefined=undef;

  if ($M2_utr_left_start){ # UTR OVERLAP UTR. As we know M1_mrna has utr in his right (otherwise "return") if we have UTR in M2 left, it is sure they overlap beetween them
    my $length_utr_M1 =  $M1_utr_right_end - $M1_utr_right_start +1;
    my $length_utr_M2 =  $M2_utr_left_end - $M2_utr_left_start + 1;
    my $dist_between_M1_and_M2 =  $M2_cds_start - $M1_cds_end - 1; # /!\ CALCUL DIFFERENT
   
    my $separting_point = undef;
    if($length_utr_M1 > $length_utr_M2){
      $separting_point = $M2_utr_left_start;
    }
    #print "($length_utr_M1 > $dist_between_M1_and_M2 and $length_utr_M2 > $dist_between_M1_and_M2)\n" if $verbose;
    # If $dist_between_M1_and_M2 = 1 we cannot share one nucleotide between two utr. Both should have it
    if($length_utr_M1 > $dist_between_M1_and_M2 and $length_utr_M2 > $dist_between_M1_and_M2 and $dist_between_M1_and_M2 > 1) {
      $separting_point =  $M1_cds_end  + int($dist_between_M1_and_M2 / 2) + 1;
      #print "separting_point  $separting_point $M2_cds_start - $M1_cds_end - 1 =  $dist_between_M1_and_M2 ".int($dist_between_M1_and_M2 / 2)."\n" if $verbose;
    }

    if($separting_point){
      print "lets shrink the UTR M1 (separting_point=$separting_point): \n";
      foreach my $l3_type (keys %{$omniscient->{'level3'}} ){               
        if ($l3_type =~ 'utr' or $l3_type =~ 'exon'){ #lets shrink it
          if(exists_keys($omniscient,('level3', $l3_type, lc($l2_M1_id)))){
            my @new_list_of_feature3;
            foreach my $feature_l3 ( sort {$a->start <=> $b->start} @{$omniscient->{'level3'}{$l3_type}{lc($l2_M1_id)}} ){
              if($feature_l3->start() >= $separting_point){ #feature is completely in the cds we remove it
                 print "we throw the feature case1A ".$feature_l3->gff_string."\n"  if $verbose;
              }
              elsif ($feature_l3->start() <= $separting_point and $feature_l3->end() >= $separting_point){
                print "we modfify the feature case1: ".$feature_l3->gff_string."\n" if $verbose;
                $feature_l3->end($separting_point-1);
                print "after modification ".$feature_l3->gff_string."\n" if $verbose;
                push(@new_list_of_feature3,  $feature_l3);
              }
              else{ #feature not concerned we keep it
                push(@new_list_of_feature3,  $feature_l3);
              }
            }
            if (scalar @new_list_of_feature3 > 0){
              @{$omniscient->{'level3'}{$l3_type}{lc($l2_M1_id)}}=@new_list_of_feature3;
            }
            else{
              print "delete it\n";
                          delete $omniscient->{'level3'}{$l3_type}{lc($l2_M1_id)};
              #clean_from_l3($omniscient, $l3_type, $l2_feature_M1);
            }
          }
        }
      }
      $utr_redefined=1;
    }
  }

  elsif ($M1_utr_right_start <= $M2_cds_end and $M1_utr_right_end >= $M2_cds_start){ #overlap CDS
    print "shrink_right utr of model 1\n";
    foreach my $l3_type (keys %{$omniscient->{'level3'}} ){               
      if ($l3_type =~ 'utr' or $l3_type =~ 'exon'){ #lets shrink it
        if(exists_keys($omniscient,('level3', $l3_type, lc($l2_M1_id)))){
          my @new_list_of_feature3;
          foreach my $feature_l3 ( sort {$a->start <=> $b->start} @{$omniscient->{'level3'}{$l3_type}{lc($l2_M1_id)}} ){
            if($feature_l3->start() >= $M2_cds_start){ #feature is completely in the cds we remove it
               print "we throw the feature case1B ".$feature_l3->gff_string."\n"  if $verbose;
            }
            elsif ($feature_l3->start() <= $M2_cds_end and $feature_l3->end() >= $M2_cds_start){
              print "we modfify the feature case1: ".$feature_l3->gff_string."\n" if $verbose;
              $feature_l3->end($M2_cds_start-1);
              push(@new_list_of_feature3,  $feature_l3);
            }
            else{ #feature not concerned we keep it
              push(@new_list_of_feature3,  $feature_l3);
            }
          }
          if (scalar @new_list_of_feature3 > 0){
            @{$omniscient->{'level3'}{$l3_type}{lc($l2_M1_id)}}=@new_list_of_feature3;
          }
          else{
            print "delete it\n";
            clean_from_l3($omniscient, $l3_type, $l2_feature_M1);
          }
        }
      }
    }
    $utr_redefined=1;
  }
  else{ #utr left de not overlap => continue
    print "Interesting, right utr overlap but not in UTR. Probably several features where overlaping and this one where more right than the other already studied \n" if $verbose;
  }


  #######################################
  # control sanity l2 and l1 location after having modified l3 location
  if($utr_redefined){
    check_mrna_positions($l2_feature_M1, $omniscient->{'level3'}{'exon'}{lc($l2_M1_id)}, $verbose);
    # print "mrna_id2 $mrna_id2".$mrna_feature2_clean->gff_string."\n";
    # check gene feature extremities
    my $l1_M1_id = $l2_feature_M1->_tag_value('Parent');
    foreach my $tag_level1 (keys %{$omniscient->{'level1'}}){
      if(exists_keys($omniscient,('level1', $tag_level1, lc($l1_M1_id)))){
        check_level1_positions($omniscient, $omniscient->{'level1'}{$tag_level1}{lc($l1_M1_id)}, $verbose);
      }
    }
  }
}



############################################################
#check if UTR left of model right is overlaping, and fix it.
sub _control_utr_from_model_right{
  my ($l2_feature_M1, $M1_cds_start, $M1_cds_end, $l2_feature_M2, $M2_cds_start, $M2_cds_end, $omniscient, $verbose)=@_;

  print "_control_utr_left_from_model_right\n";

  ####################################
  #peculiar case CDS1 and CD2 overlap. We dont touch that case
  if ($M1_cds_start <= $M2_cds_end and $M1_cds_end >= $M2_cds_start){
    print "peculiar case CDS1 and CD2 overlap. We dont touch that case\n";
    return;
  }

  my $l2_M1_id  = lc($l2_feature_M1->_tag_value('ID'));
  my $l2_M2_id  = lc($l2_feature_M2->_tag_value('ID')); 

  ##############################
  # Look at utr of the M1
  my ($M1_utr_left_start, $M1_utr_left_end, $M1_utr_right_start, $M1_utr_right_end) = _get_utrs_extremities($l2_M1_id, $M1_cds_end, $omniscient);

  ###########################
  # Look at utr of the M2
  my ($M2_utr_left_start, $M2_utr_left_end, $M2_utr_right_start, $M2_utr_right_end) = _get_utrs_extremities($l2_M2_id, $M2_cds_end, $omniscient);
  if(! $M2_utr_left_start){return;} #No UTR at left of model2, nothing to modify

  ##############################
  # CHECK RIGHT UTR OF MODEL 1 #
  ##############################
  my $utr_redefined=undef;

  ###############analysis cases
  if ($M1_utr_right_start){ # UTR OVERLAP UTR. As we know M1_mrna has utr in his right (otherwise "return") if we have UTR in M2 left, it is sure they overlap beetween them
     print "case1 ?\n";
    my $length_utr_M1 =  $M1_utr_right_end - $M1_utr_right_start +1;
    my $length_utr_M2 =  $M2_utr_left_end - $M2_utr_left_start + 1;
    my $dist_between_M1_and_M2 =  $M2_cds_start - $M1_cds_end - 1; # /!\ CALCUL DIFFERENT

    my $separting_point = undef;
    if($length_utr_M2 > $length_utr_M1){
      $separting_point = $M1_utr_right_end;
    }
    #print "($length_utr_M1 > $dist_between_M1_and_M2 and $length_utr_M2 > $dist_between_M1_and_M2)\n" if $verbose;
     # If $dist_between_M1_and_M2 = 1 we cannot share one nucleotide between two utr. Both should have it
    if($length_utr_M1 > $dist_between_M1_and_M2 and $length_utr_M2 > $dist_between_M1_and_M2 and $dist_between_M1_and_M2 > 1) {
      $separting_point =  $M1_cds_end  + int($dist_between_M1_and_M2 / 2) + 1;
      #print "separting_point $separting_point\n" if $verbose;
    }

    if($separting_point){
      print "lets shrink the UTR M1: \n";
      foreach my $l3_type (keys %{$omniscient->{'level3'}} ){               
        if ($l3_type =~ 'utr' or $l3_type =~ 'exon'){ #lets shrink it
          if(exists_keys($omniscient,('level3', $l3_type, lc($l2_M2_id)))){
            my @new_list_of_feature3;
            foreach my $feature_l3 ( sort {$a->start <=> $b->start} @{$omniscient->{'level3'}{$l3_type}{lc($l2_M2_id)}} ){
              if($feature_l3->end() <= $separting_point){ 
                 print "we throw the feature caseX1 ".$feature_l3->gff_string."\n"  if $verbose;
              }
              elsif ($feature_l3->start() <= $separting_point and $feature_l3->end() >= $separting_point){
                print "we modfify the feature caseX1: ".$feature_l3->gff_string."\n" if $verbose;
                $feature_l3->start($separting_point+1);
                push(@new_list_of_feature3,  $feature_l3);
              }
              else{ #feature not concerned we keep it
                push(@new_list_of_feature3,  $feature_l3);
              }
            }
            if (scalar @new_list_of_feature3 > 0){
              @{$omniscient->{'level3'}{$l3_type}{lc($l2_M2_id)}}=@new_list_of_feature3;
            }
            else{
              print "delete it\n";
              clean_from_l3($omniscient, $l3_type, $l2_feature_M2);
            }
          }
        }
      }
      $utr_redefined=1;
    }
  }
  elsif ($M2_utr_left_start <= $M1_cds_end and $M2_utr_left_end >= $M1_cds_start){ #overlap CDS
      print "case2 ?\n";
    #shrink_left_utr of model 2
    foreach my $l3_type (keys %{$omniscient->{'level3'}} ){               
      if ($l3_type =~ 'utr' or $l3_type =~ 'exon'){ #lets shrink it
        if(exists_keys($omniscient,('level3', $l3_type, lc($l2_M2_id)))){
          my @new_list_of_feature3;
          foreach my $feature_l3 ( sort {$a->start <=> $b->start} @{$omniscient->{'level3'}{$l3_type}{lc($l2_M2_id)}} ){
            print "$l3_type ".$feature_l3->start()."<= $M1_cds_end and ".$feature_l3->end()." >= $M1_cds_end \n";
            if($feature_l3->end() <= $M1_cds_end){ 
               print "we throw the feature case1C ".$feature_l3->gff_string."\n"  if $verbose;
            }
            elsif ($feature_l3->start() <= $M1_cds_end and $feature_l3->end() >= $M1_cds_end){
              print "we modfify the feature case1: ".$feature_l3->gff_string."\n" if $verbose;
              $feature_l3->start($M1_cds_end+1);
              push(@new_list_of_feature3,  $feature_l3);
            }
            else{ #feature not concerned we keep it
              push(@new_list_of_feature3,  $feature_l3);
            }
          }
          if (scalar @new_list_of_feature3 > 0){
            @{$omniscient->{'level3'}{$l3_type}{lc($l2_M2_id)}}=@new_list_of_feature3;
          }
          else{
            print "delete it\n";
            clean_from_l3($omniscient, $l3_type, $l2_feature_M2);
          }
        }
      }
    }
    $utr_redefined=1;
  }
  else{ #utr left de not overlap => continue
    print "Interesting2, right utr overlap but not in UTR. Probably several features where overlaping and this one where more right than the other already studied \n" if $verbose;
  }

  #######################################
  # control sanity l2 and l1 location after having modified l3 location
  if($utr_redefined){
    check_mrna_positions($l2_feature_M2, $omniscient->{'level3'}{'exon'}{lc($l2_M2_id)}, $verbose);
    # print "mrna_id2 $mrna_id2".$mrna_feature2_clean->gff_string."\n";
    # check gene feature extremities
    my $l1_M2_id = $l2_feature_M2->_tag_value('Parent');
    foreach my $tag_level1 (keys %{$omniscient->{'level1'}}){
      if(exists_keys($omniscient,('level1', $tag_level1, lc($l1_M2_id)))){
        check_level1_positions($omniscient, $omniscient->{'level1'}{$tag_level1}{lc($l1_M2_id)}, $verbose);
      }
    }
  }
}

sub _get_cds_location{
  my ($l2_feature, $omniscient)=@_;

  my $cds_start = undef;
  my $cds_end = undef;
  my $l2_id = $l2_feature->_tag_value('ID'); 
  if(exists_keys($omniscient,('level3', 'cds', lc($l2_id)))){
    sort {$a->start <=> $b->start} @{$omniscient->{'level3'}{'cds'}{lc($l2_id)}};
    $cds_start  = $omniscient->{'level3'}{'cds'}{lc($l2_id)}[0]->start; #first element of the array
    $cds_end = $omniscient->{'level3'}{'cds'}{lc($l2_id)}[$#{$omniscient->{'level3'}{'cds'}{lc($l2_id)}}]->end; #last element of the array  
  }
  else{
    print "unextpected the feature".$l2_feature->gff_string." doesnt have cds while it has a utr...!\n"; exit;
  }
  return $cds_start, $cds_end;
}

sub _get_utrs_extremities{
  my ($l2_id, $cds_end, $omniscient)=@_;
  my $verbose = undef;
  my $utr_left_start=undef;
  my $utr_left_end=undef;
  my $utr_right_start=undef;
  my $utr_right_end=undef;

  foreach my $l3_type (keys %{$omniscient->{'level3'}} ){
    print "$l3_type\n" if $verbose;
    if ($l3_type =~ 'utr'){

      if(exists_keys($omniscient,('level3', $l3_type, lc($l2_id)))){
        print "exists_keys $l3_type for $l2_id\n" if $verbose;
        sort {$a->start <=> $b->start} @{$omniscient->{'level3'}{$l3_type}{lc($l2_id)}};

        my $utr_start  = $omniscient->{'level3'}{$l3_type}{lc($l2_id)}[0]->start; #first element of the array 
        my $utr_end = $omniscient->{'level3'}{$l3_type}{lc($l2_id)}[$#{$omniscient->{'level3'}{$l3_type}{lc($l2_id)}}]->end; #last element of the array
        
         #check with utr is that 5' or 3'
        if($utr_start > $cds_end ){
          $utr_right_start = $utr_start ;
          $utr_right_end   = $utr_end;
        }
        else{
          $utr_left_start = $utr_start ;
          $utr_left_end   = $utr_end;
        }
      }
    }
  }
  return $utr_left_start, $utr_left_end, $utr_right_start, $utr_right_end;
}

sub clean_from_l3{
  my ($omniscient, $l3_type, $l2_feature) = @_;

  my $l2_id = lc($l2_feature->_tag_value('ID'));

  #Clean l3
  delete $omniscient->{'level3'}{$l3_type}{l2_id};

  #clean l2
  foreach my $l3_type (keys %{$omniscient->{'level3'}} ){   
    if(exists_keys($omniscient,('level3', $l3_type, $l2_id))){
      return;
    }
  }
  my $l2_parent = lc($l2_feature->_tag_value('Parent'));
  foreach my $l2_type (keys %{$omniscient->{'level2'}} ){
    if(exists_keys($omniscient,('level2', $l2_type, $l2_parent))){
      print Dumper($omniscient->{'level2'}{$l2_type}{$l2_parent});
      foreach my $feature (@{$omniscient->{'level2'}{$l2_type}{$l2_parent}}){
        my $id = lc($feature->_tag_value('ID'));
        if ($id eq $l2_id){
          print "removing $l2_id\n";
          delete $omniscient->{'level2'}{$l2_type}{$l2_parent};
          #remove l1
          if ($#{$omniscient->{'level2'}{$l2_type}{$l2_parent}} == 0){
            foreach my $tag (keys %{$omniscient->{'level1'}}){
              if(exists_keys($omniscient,('level1', $tag, $l2_parent))){
                delete $omniscient->{'level1'}{$tag}{$l2_parent};
              }
            }
          }
        }
      }
    }
  }
}

sub _get_right_utrs{
  my ($hash_omniscient, $l2_feature, $cds_end) = @_;

  my $l2_id = lc($l2_feature->_tag_value('ID'));

  foreach my $tag_l3 (keys %{$omniscient->{'level3'}}){
    if($tag_l3 =~ "utr"){
      if (exists_keys ($omniscient, ('level2', $tag_l3, $l2_id) ) ){          
          my @sorted_utr = sort {$a->start <=> $b->start} @{$omniscient->{'level3'}{$tag_l3}{$l2_id}};
          if ($sorted_utr[0]->start > $cds_end)
            return \@sorted_utr;
        }
      }
    }
  }
  return undef;
}

sub _get_left_utrs{
  my ($hash_omniscient, $l2_feature, $cds_start) = @_;

  my $l2_id = lc($l2_feature->_tag_value('ID'));

  foreach my $tag_l3 (keys %{$omniscient->{'level3'}}){
    if($tag_l3 =~ "utr"){
      if (exists_keys ($omniscient, ('level2', $tag_l3, $l2_id) ) ){          
          my @sorted_utr = sort {$a->start <=> $b->start} @{$omniscient->{'level3'}{$tag_l3}{$l2_id}};
          if ($sorted_utr[$#sorted_utr]->end < $cds_start)
            return \@sorted_utr;
        }
      }
    }
  }
  return undef;
}

__END__

=head1 NAME
 
gff3_sp_clip_UTRs.pl - 
This script focuses on UTR and it aims at cleaning overpredicted UTRs (e.g when annotating with RNAseq unstranded in Fungi). It will clip the left/right UTRs to avoid overlaps with other UTR/cds.
The only case where a UTR overlaping with something is not clipped, it is when the CDS of the reference mRNA is overlaping the CDS of the neighbor mRNA investigated.


=head1 SYNOPSIS

    ./gff3_sp_clip_UTRs.pl --gff annotation.gff --out=outFile 
    ./gff3_sp_clip_UTRs.pl --help

=head1 OPTIONS

=over 8

=item B<--gff> or  B<-g>

Input GFF3 file(s) used as M1erence. 


=item  B<--out>, B<--output>, B<--outfile> or B<-o>

Output gff3 containing the M1erence annotation with all the non-overlapping newly added genes from addfiles.gff.

=item B<--help> or B<-h>

Display this helpful text.

=back

=cut
