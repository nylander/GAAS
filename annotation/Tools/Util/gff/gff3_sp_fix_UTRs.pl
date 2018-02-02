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
my ($omniscient, $mRNAGeneLink) = BILS::Handler::GXFhandler->slurp_gff3_file_JD($gff);
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
my $sortBySeq = _sort_by_seq_no_strand($omniscient);
my %alreadyChecked;

foreach my $locusID ( keys %{$sortBySeq}){ # tag_l1 = gene or repeat etc...
  if ( exists_keys( $sortBySeq, ( $locusID, 'level1') ) ){
    foreach my $tag_l1 ( keys %{$sortBySeq->{$locusID}{'level1'}} ) { 

      my $to_check = clone($sortBySeq->{$locusID}{'level1'}{$tag_l1});

      # Go through location from left to right ### !!
      foreach my $id_l1 ( sort {$sortBySeq->{$locusID}{'level1'}{$tag_l1}{$a}[1] <=> $sortBySeq->{$locusID}{'level1'}{$tag_l1}{$b}[1] } keys %{$sortBySeq->{$locusID}{'level1'}{$tag_l1}} ) {

        if(! exists($alreadyChecked{$id_l1})){

          #remove itself 
          delete $to_check->{$id_l1};

          my $location = $sortBySeq->{$locusID}{'level1'}{$tag_l1}{$id_l1}; # This location will be updated on the fly

          # Go through location from left to right ### !!
          foreach my $id2_l1 ( sort {$to_check->{$a}[1] <=> $to_check->{$b}[1] } keys %{$to_check} ) {

            my $location_to_check = $to_check->{$id2_l1};

            #If location_to_check start if over the end of the M1erence location, we stop
            if($location_to_check->[1] > $sortBySeq->{$locusID}{'level1'}{$tag_l1}{$id_l1}[2]) {last;} 

            # Let's check if neighbor gene overlap at Gene LEVEL !
            if (($location->[1] <= $location_to_check->[2]) and ($location->[2] >= $location_to_check->[1])){
                ### OVERLAP ###
                #let's check at UTR level
               _check_gene_overlap_at_UTR($omniscient , $id_l1, $id2_l1, $verbose); #If contains UTR

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

# Check the presence of UTR of both sides of each gene
sub _check_gene_overlap_at_UTR{
  my  ($hash_omniscient, $gene_id, $gene_id2, $verbose)=@_;

  my $verbose=1;
  my $utr_mrna1 = undef;
  my $utr_mrna2 = undef;
  print "lets go: $gene_id, $gene_id2 \n";

  foreach my $l2_type (keys %{$hash_omniscient->{'level2'}} ){

    #check full CDS for each mRNA
    if(exists_keys($hash_omniscient,('level2', $l2_type, lc($gene_id)))){
      foreach my $mrna_feature (@{$hash_omniscient->{'level2'}{$l2_type}{lc($gene_id)}}){       

        if(exists_keys($hash_omniscient,('level2', $l2_type, lc($gene_id2)))){
            
          foreach my $mrna_feature2 (@{$hash_omniscient->{'level2'}{$l2_type}{lc($gene_id2)}}){ # from here bothe feature level2 are the same type

            #check utr if utr present
            foreach my $l3_type (keys %{$hash_omniscient->{'level3'}} ){
              if ($l3_type =~ 'utr'){
                  if(exists_keys($hash_omniscient,('level3', $l3_type, lc($mrna_feature->_tag_value('ID'))))){
                  print "There is UTR in M1\n";

                  my ($M1_cds_start, $M1_cds_end) = _get_cds_location($mrna_feature, $omniscient);
                  my ($M2_cds_start, $M2_cds_end) = _get_cds_location($mrna_feature2, $omniscient);

                  if($M2_cds_start < $M1_cds_start){ #cds of M2 is left of cds of M1, we have to flip the case
                    print "flip test1 = true\n";
                     _control_utr_from_model_right($mrna_feature2, $M2_cds_start, $M2_cds_end, $mrna_feature, $M1_cds_start, $M1_cds_end, $hash_omniscient, $verbose);
                     last;
                  }
                  else{
                    _control_utr_from_model_left($mrna_feature, $M1_cds_start, $M1_cds_end, $mrna_feature2, $M2_cds_start, $M2_cds_end, $hash_omniscient, $verbose);
                    last;
                  }
                }
              }
            }

            foreach my $l3_type (keys %{$hash_omniscient->{'level3'}} ){
              if ($l3_type =~ 'utr'){
                if(exists_keys($hash_omniscient,('level3', $l3_type, lc($mrna_feature2->_tag_value('ID'))))){
                  print "There is UTR in M2\n";

                  my ($M1_cds_start, $M1_cds_end) = _get_cds_location($mrna_feature, $omniscient);
                  my ($M2_cds_start, $M2_cds_end) = _get_cds_location($mrna_feature2, $omniscient);

                  if($M2_cds_start < $M1_cds_start){ #cds of M2 is left of cds of M1, we have to flip the case
                    print "flip test2 = true\n";
                    _control_utr_from_model_left($mrna_feature2, $M2_cds_start, $M2_cds_end, $mrna_feature, $M1_cds_start, $M1_cds_end, $hash_omniscient, $verbose);        
                     last;
                  }
                  else{
                    _control_utr_from_model_right($mrna_feature, $M1_cds_start, $M1_cds_end, $mrna_feature2, $M2_cds_start, $M2_cds_end, $hash_omniscient, $verbose);
                    last;
                  }
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
      print "lets shrink the UTR M1: \n";
      foreach my $l3_type (keys %{$omniscient->{'level3'}} ){               
        if ($l3_type =~ 'utr' or $l3_type =~ 'exon'){ #lets shrink it
          if(exists_keys($omniscient,('level3', $l3_type, lc($l2_M1_id)))){
            my @new_list_of_feature3;
            foreach my $feature_l3 ( sort {$a->start <=> $b->start} @{$omniscient->{'level3'}{$l3_type}{lc($l2_M1_id)}} ){
              if($feature_l3->start() >= $separting_point){ #feature is completely in the cds we remove it
                 print "we throw the feature case1 ".$feature_l3->gff_string."\n"  if $verbose;
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
            @{$omniscient->{'level3'}{$l3_type}{lc($l2_M1_id)}}=@new_list_of_feature3;
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
               print "we throw the feature case1 ".$feature_l3->gff_string."\n"  if $verbose;
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
          @{$omniscient->{'level3'}{$l3_type}{lc($l2_M1_id)}}=@new_list_of_feature3;
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
  # Look at utr of the M1erence
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
            @{$omniscient->{'level3'}{$l3_type}{lc($l2_M2_id)}}=@new_list_of_feature3;
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
            print "kkk\n";
            print "$l3_type ".$feature_l3->start()."<= $M1_cds_end and ".$feature_l3->end()." >= $M1_cds_end \n";
            if($feature_l3->end() <= $M1_cds_end){ 
               print "we throw the feature case1 ".$feature_l3->gff_string."\n"  if $verbose;
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
          @{$omniscient->{'level3'}{$l3_type}{lc($l2_M2_id)}}=@new_list_of_feature3;
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

      if(exists_keys($omniscient,('level3', $l3_type, $l2_id))){
        print "exists_keys $l3_type for $l2_id\n" if $verbose;
        sort {$a->start <=> $b->start} @{$omniscient->{'level3'}{$l3_type}{lc($l2_id)}};
        #print Dumper($omniscient->{'level3'}{$l3_type}{lc($l2_id)});
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

#
# Sort by locusID !!!!
# L1 => LocusID->level->typeFeature->ID =[ID,start,end]
# L2 and L3 => LocusID->level->typeFeature->Parent->ID = [ID,start,end]
#
#
sub _sort_by_seq_no_strand{
  my ($omniscient) = @_;

  my %hash_sortBySeq;
  
    foreach my $tag_level1 (keys %{$omniscient->{'level1'}}){
      foreach my $level1_id (keys %{$omniscient->{'level1'}{$tag_level1}}){
        my $level1_feature = $omniscient->{'level1'}{$tag_level1}{$level1_id};
        my $ID = $level1_feature->_tag_value('ID');
        my $position_l1=$level1_feature->seq_id;

        $hash_sortBySeq{$position_l1}{"level1"}{$tag_level1}{$level1_id} = [$ID, int($level1_feature->start), int($level1_feature->end)];
        }
  }
  return \%hash_sortBySeq;
}

__END__

=head1 NAME
 
gff3_sp_fix_UTRs.pl - 
This script focuses on UTR and it aims at cleaning overpredicted UTRs (e.g when annotating with RNAseq unstranded in Fungi). It will clip the left/right UTRs to avoid overlaps with other UTR/cds.
The only case where a UTR overlaping with something is not clipped, it is when the CDS of the reference mRNA is overlaping the CDS of the neighbor mRNA investigated.


=head1 SYNOPSIS

    ./gff3_sp_fix_UTRs.pl --gff annotation.gff --out=outFile 
    ./gff3_sp_fix_UTRs.pl --help

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
