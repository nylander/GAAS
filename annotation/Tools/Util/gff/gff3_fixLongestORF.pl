#!/usr/bin/perl

## if IUPAC:
## We consider a stop only if we are sure it is one
## CDS can contains putative stop codon (but not sure stop one like YAA that can be TAA or CAA).
## We consider a start even if is not sure like AYG that can be ATG or ACG

use Carp;
use Clone 'clone';
use strict;
use File::Basename;
use Getopt::Long;
use Statistics::R;
use Pod::Usage;
use Data::Dumper;
use List::MoreUtils qw(uniq);
use Bio::Tools::GFF;
use Bio::DB::Fasta;
#use Bio::Seq;
use Bio::SeqIO;
use BILS::Handler::GXFhandler qw(:Ok);
use BILS::Handler::GFF3handler qw(:Ok);
use BILS::Plot::R qw(:Ok);

my $SIZE_OPT=21;

my $header = qq{
########################################################
# BILS 2015 - Sweden                                   #  
# jacques.dainat\@bils.se                               #
# Please cite BILS (www.bils.se) when using this tool. #
########################################################
};

my $outfile = undef;
my $gff = undef;
my $model_to_test = undef;
my $file_fasta=undef;
my $split_opt=undef;
my $help= 0;

if ( !GetOptions(
    "help|h" => \$help,
    "gff=s" => \$gff,
    "fasta|fa=s" => \$file_fasta,
    "split|s" => \$split_opt,
    "m|model=s" => \$model_to_test,
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
 
if ( ! (defined($gff)) or !(defined($file_fasta)) ){
    pod2usage( {
           -message => "$header\nAt least 2 parameter is mandatory:\nInput reference gff file (--gff) and Input fasta file (--fasta)\n\n",
           -verbose => 0,
           -exitval => 1 } );
}

######################
# Manage output file #
my $gffout;
my $gffout2;
my $gffout3;
#my $gffout4;
if ($outfile) {
  $outfile=~ s/.gff//g;
open(my $fh, '>', $outfile."-intact.gff") or die "Could not open file '$outfile' $!";
  $gffout= Bio::Tools::GFF->new(-fh => $fh, -gff_version => 3 );
open(my $fh2, '>', $outfile."-only_modified.gff") or die "Could not open file '$outfile' $!";
  $gffout2= Bio::Tools::GFF->new(-fh => $fh2, -gff_version => 3 );
open($gffout3, '>', $outfile."-report.txt") or die "Could not open file '$outfile' $!";
#open(my $fh3, '>', $outfile."-pseudogenes.gff") or die "Could not open file '$outfile' $!";
#  $gffout4= Bio::Tools::GFF->new(-fh => $fh3, -gff_version => 3 );
}
else{
  $gffout = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 3);
}

my %ListModel;
if(!($model_to_test)){
  $ListModel{1}=0;
  $ListModel{2}=0;
  $ListModel{3}=0;
  $ListModel{4}=0;
}else{
  my @fields= split(',', $model_to_test);
  foreach my $field (@fields){
    if($field =~ m/^[01234]$/){
      $ListModel{$field}=0;
    }else{
      print "This model $field is not known. Must be an Integer !\n";exit;      
    }
  }
}

                #####################
                #     MAIN          #
                #####################


######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) = BILS::Handler::GFF3handler->slurp_gff3_file_JD($gff);
print ("GFF3 file parsed\n");


####################
# index the genome #
my $db = Bio::DB::Fasta->new($file_fasta);
print ("Genome fasta parsed\n");

####################
my $pseudo_threshold=70;
#counters
my $counter_case21=0;
my $geneCounter=0;
my $mRNACounter=0;
my $mRNACounter_fixed=0;
#my $mrna_pseudo_suspected=0;
#my $gene_pseudo_suspected=0;
#my $mrna_pseudo_removed=0;
#my $gene_pseudo_removed=0;
my $special_or_partial_mRNA=0;

my %omniscient_modified_gene;
#my %omniscient_pseudogene;
my @modified_gene_list;
my @intact_gene_list;

foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
  foreach my $gene_id_tag_key (keys %{$hash_omniscient->{'level1'}{$primary_tag_key_level1}}){
      my $gene_feature=$hash_omniscient->{'level1'}{$primary_tag_key_level1}{$gene_id_tag_key};

    my $one_ORFmodified="no";
    #my $mrna_pseudo=0;
    #my @list_mrna_pseudo;
    my $one_level2_modified; # check if one of the level2 feature will be modified
    my $number_mrna=0;

    foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
      foreach my $level2_feature ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$gene_id_tag_key}}) {
       
        my $ORFmodified="no";
        $number_mrna=$#{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$gene_id_tag_key}}+1;
       
        # get level2 id
        my @values = $level2_feature->get_tag_values('ID');       
        my $id_level2 = lc(shift @values) ;
        
        ##############################
        #If it's a mRNA = have CDS. #
        if ( exists ($hash_omniscient->{'level3'}{'cds'}{$id_level2} ) ){
                  
          ##############
          # Manage CDS #
          my @cds_feature_list = sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{'cds'}{$id_level2}}; # be sure that list is sorted
          my ($cdsExtremStart, $cds_dna_seq, $cdsExtremEnd) = concatenate_feature_list(\@cds_feature_list);
          #create the cds object
          my $cds_obj = Bio::Seq->new(-seq => $cds_dna_seq, -alphabet => 'dna' );         
          #Reverse the object depending on strand
          if ($level2_feature->strand == -1 or $level2_feature->strand eq "-"){
            $cds_obj = $cds_obj->revcom();
          }
          #translate cds in protein
          my $original_prot_obj = $cds_obj->translate() ; #codontable_id by default=1 (Vertebrates). IUPAC => STOP codon even if not sure ...
          my $cds_prot=$original_prot_obj->seq;
          my $originalProt_size=length($cds_prot);

          ################################################
          # mRNA: extract the concatenated exon sequence #
          my @exons_features = sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{'exon'}{$id_level2}};
          my ($exonExtremStart, $mrna_seq, $exonExtremEnd) = concatenate_feature_list(\@exons_features);
          #create the mrna object
          my $mrna_obj = Bio::Seq->new(-seq => $mrna_seq, -alphabet => 'dna' );
  #        print $mrna_obj->seq."\n"; 
          #Reverse complement according to strand
          if ($level2_feature->strand == -1 or $level2_feature->strand eq "-"){
            $mrna_obj = $mrna_obj->revcom();
          }
          
          #######################
          # Get the longest ORF ## record ORF = start, end (half-open), length, and frame
          my $longest_ORF_prot_obj;
          my $orf_cds_region;
          my ($longest_ORF_prot_objM, $orf_cds_regionM) = translate_JD($mrna_obj, 
                                                                      -nostartbyaa => 'L',
                                                                      -orf => 'longest');
          my ($longest_ORF_prot_objL, $orf_cds_regionL) = translate_JD($mrna_obj, 
                                                                      -nostartbyaa => 'M',
                                                                      -orf => 'longest');
          if($longest_ORF_prot_objL->length()+$SIZE_OPT > $longest_ORF_prot_objM ){ # In a randomly generated DNA sequence with an equal percentage of each nucleotide, a stop-codon would be expected once every 21 codons. Deonier et al. 2005
            $longest_ORF_prot_obj=$longest_ORF_prot_objL;                    # As Leucine L (9/100) occur more often than Metionine M (2.4) JD arbitrary choose to use the L only if we strength the sequence more than 21 AA. Otherwise we use M start codon.
            $orf_cds_region=$orf_cds_regionL;
            $counter_case21++;
          }
          else{

            $longest_ORF_prot_obj=$longest_ORF_prot_objM;
            $orf_cds_region=$orf_cds_regionM;

          }
    #       print Dumper($orf_cds_region)."\n";
          # set real start and stop to orf
          my $realORFstart;
          my $realORFend;
          # change the start for negative strand
          if ($level2_feature->strand == -1 or $level2_feature->strand eq "-"){
            $orf_cds_region->[0]=(length($mrna_seq) - $orf_cds_region->[1]); 
          }
          #calcul the real start end stop of cds in genome
   #       print Dumper($orf_cds_region)."\n".$mrna_obj->seq."\n";  
          ($realORFstart, $realORFend) = calcul_real_orf_end_and_start($orf_cds_region, \@exons_features);
  #        print "$id_level2 $realORFstart $realORFend\n";         
          #save the real start and stop
          $orf_cds_region->[0]=$realORFstart;
          $orf_cds_region->[1]=$realORFend;       

#############
# Tests     #
#############

          ########################
          # prediction is longer #
          if($longest_ORF_prot_obj->length() > $originalProt_size){

#Model1     ###############################################
            # sequence original is part of new prediction #
            if (index($longest_ORF_prot_obj->seq,$cds_prot) != -1){
              if ( exists($ListModel{1}) ){ 
                if(!(($longest_ORF_prot_obj->seq =~ m/^X/) and ($longest_ORF_prot_obj->length() < $originalProt_size+$SIZE_OPT))){ #avoid case of ambigous methionine (Written X) -> Need to be over 21 AA to decide ok is longer and can be a M

                  print "mrNA $id_level2 => model1. prot original is shorter:\noriginal:$cds_prot\ncdsStart $cdsExtremStart - cdsEnd $cdsExtremEnd\n".
                  $longest_ORF_prot_obj->seq."\n\n";
                  $ListModel{1}++;
                  modify_gene_model($hash_omniscient, \%omniscient_modified_gene, $gene_feature, $gene_id_tag_key, $level2_feature, $id_level2, \@exons_features, \@cds_feature_list, $cdsExtremStart, $cdsExtremEnd, $realORFstart, $realORFend, 'model1', $gffout);
                  $ORFmodified="yes";

                }
              }
            }
            #################################################
            # protein original and predicted are different
            else{

              #########################################
#Model2       # Prediction don't overlap original CDS #
              if( ($realORFend < $cdsExtremStart) or ($realORFstart > $cdsExtremEnd) ){ 
                my $model;
  
                if( exists($ListModel{2}) ){
                   print "mrNA $id_level2 $gene_id_tag_key => model1\n";
                  $ListModel{2}++;
                  $model=1;
                  if($split_opt){
                    split_gene_model(\@intact_gene_list, $hash_omniscient, \%omniscient_modified_gene, $gene_feature, $gene_id_tag_key, $level2_feature, $id_level2, \@exons_features, \@cds_feature_list, $cdsExtremStart, $cdsExtremEnd, $realORFstart, $realORFend, 'model2', $gffout);
                  }
                  else{
                    modify_gene_model($hash_omniscient, \%omniscient_modified_gene, $gene_feature, $gene_id_tag_key, $level2_feature, $id_level2, \@exons_features, \@cds_feature_list, $cdsExtremStart, $cdsExtremEnd, $realORFstart, $realORFend, 'model2', $gffout);
                  }
                  $ORFmodified="yes";
                }
              } # End they don't overlap
             
              #########################################
              # Prediction Overlap original CDS #
              else { # They overlap
#Model3         ###############
                # original protein and predicted one are different; the predicted one is longest, they overlap each other.
                if( exists($ListModel{3}) ){ 
                  $ListModel{3}++;
                  print "mrNA $id_level2 $gene_id_tag_key => model3\n";
                  
                  modify_gene_model($hash_omniscient, \%omniscient_modified_gene, $gene_feature, $gene_id_tag_key, $level2_feature, $id_level2, \@exons_features, \@cds_feature_list, $cdsExtremStart, $cdsExtremEnd, $realORFstart, $realORFend, 'model3', $gffout);
                  $ORFmodified="yes";
                }
              }
            } 
          }# End prediction longer

          ###########################
          # The real ORF looks to be shorter than the one originaly described ! Selenocysteine ? pseudogene ? Or just case where prediction begin by L instead of M (correct !) or begin by XXXXXX
          elsif($longest_ORF_prot_obj->length() < $originalProt_size){

#Model4     ###############
            # /!\ Here we compare the CDS traduction (traduct in IUPAC) against longest CDS in mRNA IUPAC modified to take in account for stops codon only those that are trustable only (TGA, TAR...).
            if( exists($ListModel{4}) ){ 
              $ListModel{4}++;
              print "mrNA $id_level2 $gene_id_tag_key => model4\n";
              #print "Original: ".$original_prot_obj->seq."\n";
              #print "longestl: ".$longest_ORF_prot_obj->seq."\n";
              # contains stop codon but not at the last position
              if( (index($original_prot_obj->seq, '*') != -1 ) and (index($original_prot_obj->seq, '*') != length($original_prot_obj->seq)-1) ){
           ## Pseudogene THRESHOLD ##              
  #              my $threshold_size=(length($original_prot_obj->seq)*$pseudo_threshold)/100; #70% of the original size
  #              if(length($longest_ORF_prot_obj->seq) <  $threshold_size){ # inferior to threshold choosen, we suspect it to be a pseudogene
  #                print Dumper($original_prot_obj);
  #                print Dumper($longest_ORF_prot_obj);
  #                $mrna_pseudo++;
  #                push(@list_mrna_pseudo, $id_level2);
  #              }
  #              else{ #remodelate a shorter gene
                  modify_gene_model($hash_omniscient, \%omniscient_modified_gene, $gene_feature, $gene_id_tag_key, $level2_feature, $id_level2, \@exons_features, \@cds_feature_list, $cdsExtremStart, $cdsExtremEnd, $realORFstart, $realORFend, 'model4', $gffout);
                  $ORFmodified="yes";
  #              }
              }# Doesn't contain stop in the middle of the sequence
              else{$special_or_partial_mRNA++;}
            }
          }


        } # End there is a CDS
        if($ORFmodified eq "yes"){
          $one_ORFmodified="yes";
          $mRNACounter_fixed++; # Count only mRNA modified
        }
      } # End foreach mRNA

  #    if($mrna_pseudo > 0){
        # all mRNA are pseudogene, we change the gene status to pseudogenes.
  #      if($mrna_pseudo == $number_mrna){
  #        $mrna_pseudo_suspected=$mrna_pseudo_suspected+$number_mrna;
  #        $gene_pseudo_suspected++;
  #       $gene_feature->primary_tag('pseudogene'); 
          #transfert the gene and sub-feature to the omniscient_pseudogene hash
  #        my @level1_list=($gene_id_tag_key);
  #        fill_omniscient_from_other_omniscient_level1_id(\@level1_list, $hash_omniscient, \%omniscient_pseudogene); # If already exists in omniscient_modified_gene, it will be replaced by the modified one
  #      }
        #only some of the isoform are pseudo... we remove them
  #      else{
  #        $mrna_pseudo_removed=$mrna_pseudo_removed+$mrna_pseudo;
  #        $gene_pseudo_removed++;
  #        my @tag_list=('all');
  #        my @id_list=($gene_id_tag_key);
  #        remove_element_from_omniscient(\@id_list, \@list_mrna_pseudo, $hash_omniscient, 'level2', 'false', \@tag_list);
  #        remove_tuple_from_omniscient(\@list_mrna_pseudo, $hash_omniscient, 'level3', 'false', \@tag_list);
  #        remove_element_from_omniscient(\@id_list, \@list_mrna_pseudo, \%omniscient_modified_gene, 'level2', 'false', \@tag_list);
  #        remove_tuple_from_omniscient(\@list_mrna_pseudo, \%omniscient_modified_gene, 'level3', 'false', \@tag_list);
  #        print "@list_mrna_pseudo has been removed because are isoform containing stop codon\n";
  #      }
  #    }

      if($one_ORFmodified eq "yes"){
        $geneCounter++;
        $mRNACounter=$mRNACounter+$number_mrna; #add all the mRNA if at least one modified
        #save remodelate gene name
        push(@modified_gene_list, $gene_id_tag_key);
      }
      else{push(@intact_gene_list, $gene_id_tag_key);}
    }
  }
}
###
# Fix frame
fil_cds_frame(\%omniscient_modified_gene);
#fil_cds_frame(\%omniscient_pseudogene);
fil_cds_frame($hash_omniscient);

########
# Print results
if ($outfile) {
  #print all in file1
  print_omniscient_from_level1_id_list($hash_omniscient, \@intact_gene_list, $gffout); #print intact gene to the file
  print_omniscient(\%omniscient_modified_gene, $gffout2); #print gene modified in file 
  #print_omniscient(\%omniscient_pseudogene, $gffout4); #print putative pseudogene in file
 
}
else{
  #print_omniscient_from_level1_id_list($hash_omniscient, \@intact_gene_list, $gffout); #print gene intact
  print_omniscient(\%omniscient_modified_gene, $gffout); #print gene modified
  #print_omniscient(\%omniscient_pseudogene, $gffout); #print putative pseudogene
}

#END
my $string_to_print="Results:\n";
$string_to_print .= "$geneCounter genes has been modified. These gene has  $mRNACounter mRNA, and among them  $mRNACounter_fixed had their ORF fixed.\n";
if (exists ($ListModel{1})){
  $string_to_print .= "$ListModel{1} model1: Prediction(s) contains the orignal prediction but is longer.\n";
}
if (exists ($ListModel{2})){
  $string_to_print .= "$ListModel{2} model2: Longest ORF found non-overlaping the original one.";
  if ($split_opt){
    $string_to_print .= " Thus, sequences have been split en two different genes (Consequently $ListModel{2} new genes has been created";
    }
     $string_to_print .= "\n";
}
if (exists ($ListModel{3})){
  $string_to_print .= "$ListModel{3} model3: sequences have been re-shaped/re-modeled (Longest ORF found overlaping the original one but doesn't contain it.)\n";
}
if (exists ($ListModel{4})){
  my $withStop=$ListModel{4}-$special_or_partial_mRNA;
  #my $withStop_butstillgene=$ListModel{4}-($mrna_pseudo_suspected+$mrna_pseudo_removed)-$special_or_partial_mRNA;
  $string_to_print .="$ListModel{4} model4: The new prediction was shorter than the original:\nAmong them, $special_or_partial_mRNA are partial (begining or finishing by NNNN or XXXX) or begining by L (It is a correct possibility), so we don't touch them.\n".
  "In other hand, $withStop are shorter due to the presence of stop codon.".#" The threshold to declare them as a pseudogene (comparing to the original size) is $pseudo_threshold percent.\n".
 # "According to this threshold, we change the gene status (prinary_tag) of $gene_pseudo_suspected genes (corresponding to $mrna_pseudo_suspected mRNA) to pseudogene.\n".
 # "According to this threshold, we suspect $gene_pseudo_suspected genes to be pseudogenes (corresponding to $mrna_pseudo_suspected mRNA). So they habe been reported in a secpific output file.\n".
 # "$withStop_butstillgene mRNA(s) containing stop but over this treshold has been re-modelate.\n".
 " They have been remodeleted.";

 # "Moreover, $mrna_pseudo_removed putative pseudo mRNA isoforms have been removed because the gene has as well non-pseudo mRNA.\n";
}

$string_to_print .="\n/!\\Remind:\n L and M are AA are possible start codons.\nParticular case: If we have a triplet as WTG, AYG, RTG, RTR or ATK it will be seen as a possible Methionine codon start (it's a X aa)\n".
"An arbitrary choisce has been done: The longer translate can begin by a L only if it's longer by 21 AA than the longer translate beginning by M. It's happened $counter_case21 times here.\n";

print $string_to_print;
if($outfile){
  print $gffout3 $string_to_print
}
print "Bye Bye.\n";
#######################################################################################################################
        ####################
         #     METHODS    #
          ################
           ##############
            ############
             ##########
              ########
               ######
                ####
                 ##
sub modify_gene_model{

  my ($hash_omniscient, $omniscient_modified_gene, $gene_feature, $gene_id_tag_key, $level2_feature, $id_level2, $exons_features, $cds_feature_list, $cdsExtremStart, $cdsExtremEnd, $realORFstart, $realORFend, $model, $gffout)=@_;

                  ###############################################
                  # modelate level3 features for new prediction #
                  my ($new_pred_utr5_list, $new_pred_cds_list, $new_pred_utr3_list) = modelate_utr_and_cds_features_from_exon_features_and_cds_start_stop($exons_features, $realORFstart, $realORFend);
                  
                  #########
                  #RE-SHAPE last/first exon if less than 3 nucleotides (1  or 2 must be romved) when the CDS finish 1 or 2 nuclotide before... because cannot be defined as UTR
                  shape_exon_extremity($exons_features, $new_pred_cds_list);
                  my ($new_pred_utr5_list, $variable_not_needed, $new_pred_utr3_list) = modelate_utr_and_cds_features_from_exon_features_and_cds_start_stop($exons_features, $realORFstart, $realORFend);

                  #############################################################
                  #  Remove ancient cds
                  my @tag_list=('exon');
                  my @id_list=($id_level2);
                  remove_tuple_from_omniscient(\@id_list, $hash_omniscient, 'level3', 'false', \@tag_list);

                  ####################
                  # Add new CDS/UTRs
                  foreach my $cds_feature (@$new_pred_cds_list){
                    push (@{$hash_omniscient->{'level3'}{'cds'}{$id_level2}}, $cds_feature);
                  }
                  foreach my $utr5_feature (@$new_pred_utr5_list){
                    push (@{$hash_omniscient->{'level3'}{'five_prime_utr'}{$id_level2}}, $utr5_feature);
                  }
                  foreach my $utr3_feature (@$new_pred_utr3_list){
                    push (@{$hash_omniscient->{'level3'}{'three_prime_utr'}{$id_level2}}, $utr3_feature);
                  }
                  $level2_feature->add_tag_value('orfix', $model);

                  check_start_end_of_mrna_feature($level2_feature, $exons_features);
                  check_start_end_of_gene_feature($hash_omniscient, $gene_id_tag_key);

                  #transfert the gene and sub-feature to the omniscient_modified_gene hash
                  my @level1_list=($gene_id_tag_key);
                  fill_omniscient_from_other_omniscient_level1_id(\@level1_list, $hash_omniscient, $omniscient_modified_gene); # If already exists in omniscient_modified_gene, it will be replaced by the modified one
              
}
############ /!\
# P.S: To be perfect, when a gene is newly created, we should verify if it is not created where another one has already been created. If yes, the should be linked together !!
############
sub split_gene_model{

   my ($intact_gene_list, $hash_omniscient, $omniscient_modified_gene, $gene_feature, $gene_id_tag_key, $level2_feature, $id_level2, $exons_features, $cds_feature_list, $cdsExtremStart, $cdsExtremEnd, $realORFstart, $realORFend, $model, $gffout)=@_;

      my $numberOfNewGene=1;

      my @values=$gene_feature->get_tag_values('ID');
      my $realGeneName=shift(@values);

                  ######################
                  # Recreate exon list #
                  my $bolean_original_is_first;
                  my $first_end;
                  my $second_start;
                  #if new prediction after on the sequence
                  if($realORFstart >= $cdsExtremEnd){
                    $bolean_original_is_first="true";
                    $first_end=$cdsExtremEnd;
                    $second_start=$realORFstart;
                  }else{ # ($realORFend < $cdsExtremStart)
                    $bolean_original_is_first="false";
                    $first_end=$realORFend;
                    $second_start=$cdsExtremStart;
                  }
                  my ($newOrignal_exon_list, $newPred_exon_list) = create_two_exon_lists($exons_features,$first_end,$second_start,$bolean_original_is_first);

        ####################################
        # Remodelate ancient gene
        ####################################

                  #############################################################
                  #  Remove all level3 feature execept cds
                  my @tag_list=('cds');
                  my @id_list=($id_level2);
                  remove_tuple_from_omniscient(\@id_list, $hash_omniscient, 'level3', 'false', \@tag_list);
                  #############
                  # Recreate original exon 
                  @{$hash_omniscient->{'level3'}{'exon'}{$id_level2}}=@$newOrignal_exon_list;

                  #########
                  #RE-SHAPE last/first exon if less than 3 nucleotides (1  or 2 must be romved) when the CDS finish 1 or 2 nuclotide before... because cannot be defined as UTR
                  shape_exon_extremity($newOrignal_exon_list,$cds_feature_list);
                 
                  ########
                  # calcul utr 
                  my ($original_utr5_list, $variable_not_needed, $original_utr3_list) = modelate_utr_and_cds_features_from_exon_features_and_cds_start_stop($newOrignal_exon_list, $cdsExtremStart, $cdsExtremEnd);

                  #########
                  #RE-SHAPE mrna extremities
                  check_start_end_of_mrna_feature($level2_feature, $newOrignal_exon_list);
                  $level2_feature->add_tag_value('orfix',$model);
                  #########
                  #RE-SHAPE gene model
                  if (must_be_a_new_gene($hash_omniscient, $gene_id_tag_key, $id_level2, $level2_feature)){
                    ## create a new gene
                    my $new_gene_id="new_".$realGeneName."-".$numberOfNewGene;
                    $numberOfNewGene++;
                    my $new_gene_feature = Bio::SeqFeature::Generic->new(-seq_id => $level2_feature->seq_id, -source_tag => $level2_feature->source_tag, -primary_tag => 'gene' , -start => $level2_feature->start,  -end => $level2_feature->end, -frame => $level2_feature->frame, -strand => $level2_feature->strand , -tag => { 'ID' => $new_gene_id }) ;
                    create_or_replace_tag($level2_feature,'Parent',$new_gene_id);

                    # append new gene in omniscient_modified_gene
                    my @level1_list=($new_gene_feature);
                    my @level2_list=($level2_feature);
                    my @level3_list=(@$newOrignal_exon_list, @$cds_feature_list, @$original_utr5_list, @$original_utr3_list);

                    append_omniscient($omniscient_modified_gene, \@level1_list, \@level2_list, \@level3_list);
                  }
                  else{ # keep the original gene model that we modified
                    ## check shape of original gene
                    check_start_end_of_gene_feature($hash_omniscient, $gene_id_tag_key);
                    $level2_feature->add_tag_value('orfix',$model);
                    # append gene modified in omniscient_modified_gene
                    my @level1_list=($gene_id_tag_key);
                    fill_omniscient_from_other_omniscient_level1_id(\@level1_list, $hash_omniscient, $omniscient_modified_gene); # If already exists in omniscient_modified_gene, it will be replaced by the modified one
                  
                  }
        ###################################
        # Remodelate New Prediction
        ###################################
                  ###############################################
                  # modelate level3 features for new prediction #
                  my ($new_pred_utr5_list, $new_pred_cds_list, $new_pred_utr3_list) = modelate_utr_and_cds_features_from_exon_features_and_cds_start_stop($newPred_exon_list, $realORFstart, $realORFend);

                  ####################################
                  #RE-SHAPE last/first exon if less than 3 nucleotides (1  or 2 must be romved) when the CDS finish 1 or 2 nuclotide before... because cannot be defined as UTR
                  shape_exon_extremity($newPred_exon_list, $new_pred_cds_list);  

                  ######################################################
                  # Modelate gene and mRNA features for new prediction #
                  my @values = $newPred_exon_list->[0]->get_tag_values('Parent');
                  my $transcript_id = shift @values;
                  my $new_mRNA_feature = Bio::SeqFeature::Generic->new(-seq_id => $newPred_exon_list->[0]->seq_id, -source_tag => $newPred_exon_list->[0]->source_tag, -primary_tag => 'mRNA' , -start => $newPred_exon_list->[0]->start,  -end => $newPred_exon_list->[$#{$newPred_exon_list}]->end, -frame => $newPred_exon_list->[0]->frame, -strand => $newPred_exon_list->[0]->strand , -tag => { 'ID' => $transcript_id , 'Parent' => $realGeneName }) ;

                  my @level1_list;
                  my @level2_list;
                  my @level3_list;

                  #$numberOfNewGene == 1 mean we already use the gene container. So in the case where we have oly one mRNA, the split will give 2 mRNA. One is linked to the original gene container (done before)
                  # The second must be linked to a new gene container. So, even if must_be_a_new_gene method say no, we must create it because the original one has been already used.
                  my $create_a_new_gene=must_be_a_new_gene($hash_omniscient, $gene_id_tag_key, $transcript_id, $new_mRNA_feature);
                  if ( ($#{$hash_omniscient->{'level2'}{'mrna'}{$gene_id_tag_key}} == 0) and $numberOfNewGene == 1){ $create_a_new_gene="true";}
                  if ( $create_a_new_gene ){ 
                    my $new_gene_id="new_".$realGeneName."-".$numberOfNewGene;
                    create_or_replace_tag($new_mRNA_feature, 'Parent', $new_gene_id);
                    my $new_gene_feature = Bio::SeqFeature::Generic->new(-seq_id => $newPred_exon_list->[0]->seq_id, -source_tag => $newPred_exon_list->[0]->source_tag, -primary_tag => 'gene' , -start => $newPred_exon_list->[0]->start,  -end => $newPred_exon_list->[$#{$newPred_exon_list}]->end, -frame => $newPred_exon_list->[0]->frame, -strand => $newPred_exon_list->[0]->strand , -tag => { 'ID' => $new_gene_id } , 'orfix' => $model) ;
                    @level1_list=($new_gene_feature);
                    @level2_list=($new_mRNA_feature);
                  }
                  else{ #the new mRNA still overlap an isoform. So we keep the link with the original gene  
                    # append new gene in omniscient_modified_gene
                    check_start_end_of_gene_feature($hash_omniscient, $gene_id_tag_key);
                    @level1_list=($gene_feature);
                    @level2_list=($new_mRNA_feature);
                  }
                  @level3_list=(@$newPred_exon_list, @$new_pred_cds_list, @$new_pred_utr5_list, @$new_pred_utr3_list);
                  append_omniscient($omniscient_modified_gene, \@level1_list, \@level2_list, \@level3_list); # If already exists , no replacement

                  if ($numberOfNewGene > 1){
                    #remove the mRNA from original omnicient (because the two mRNAs form the splited one are no linked to the original gene
                    # but are now linked to newly created gene features). The same for all linked level 3 features, so we remove them.
                    my @tag_list=('all');
                    my @id_list=($id_level2);

                    remove_tuple_from_omniscient(\@id_list, $hash_omniscient, 'level3', 'false', \@tag_list);
                    my @id_list=($gene_id_tag_key);my @id_list2=($id_level2);
                    remove_element_from_omniscient(\@id_list, \@id_list2, $hash_omniscient, 'level2', 'false', \@tag_list);
                    #reshape end and start
                    check_start_end_of_gene_feature($hash_omniscient, $gene_id_tag_key);
                    push(@{$intact_gene_list}, $gene_id_tag_key);
                  }


}

# Yes if mRNA doesnt overlap an other existing isoform
sub must_be_a_new_gene{
  my ($hash_omniscient, $gene_id, $id_level2, $level2_feature)=@_;
 
  my $result="true";
  my @list_mrna=@{$hash_omniscient->{'level2'}{'mrna'}{$gene_id}};

  if($#list_mrna > 0){ #more than only one mrna
    foreach my $mrna (@list_mrna){
      # get level2 id
      my @values = $mrna->get_tag_values('ID');       
      my $mrna_id = lc(shift @values);
      if(lc($id_level2) ne lc($mrna_id)){ # we dont check mrna against itself
        #Now check if overlap
        if( ($level2_feature->start <= $mrna->end) and ($level2_feature->end >= $mrna->start) ){ # if it overlaps
          $result=undef;last;
        }
      }
    }
  }
  else{$result=undef;} #only one mRNA

  return $result;
}

sub shape_exon_extremity{
  #exon_features is a sorted list
  #cds_features is a sorted list

  my ($exon_features,$cds_features)=@_;

   #test between first exon and first cds
   if( (abs($cds_features->[0]->start - $exon_features->[0]->start) < 3) and (abs($cds_features->[0]->start - $exon_features->[0]->start) > 0) ){ #We have to shape the exon start. We don't want a non multiple of 3 inferior to 3

      $exon_features->[0]->start($cds_features->[0]->start);
#      print "start reshaped\n";
   }
   #test between last exon and last cds
   if(abs($exon_features->[$#{ $exon_features }]->end - $cds_features->[$#{ $cds_features }]->end ) < 3){  #We have to shape the exon end
      $exon_features->[$#{ $exon_features }]->end($cds_features->[$#{ $cds_features }]->end);
#      print "end reshaped\n";
   }
}

sub calcul_real_orf_end_and_start{
  #exons_features is sorted
  my ($orf_cds_region, $exons_features)=@_;

  my $realORFstart;
  my $realORFend;

  my $orf_start=$orf_cds_region->[0]; # get start to begin
  my $orf_length=$orf_cds_region->[2]; # get start to begin

  my $first="yes";
  my $total_exon_length=0;
  my $total_exon_length_previous_round=0;
  my $mapped_length=0;
  my $mapped_length_total=0;
  my $the_rest_to_map=0; 

  foreach my $exon_feature (@$exons_features){   
    # Allows to follow the path on mRNA
    my $exon_length=($exon_feature->end - $exon_feature->start)+1;
    $total_exon_length_previous_round=$total_exon_length;
    $total_exon_length=$total_exon_length+$exon_length;
    # Allows to follow the path on the CDS
    $mapped_length_total=$mapped_length_total+$mapped_length;
    $the_rest_to_map=$orf_length-$mapped_length_total;
    # exon overlap CDS
    if($total_exon_length >= $orf_start){ #they begin to overlap
      if($first eq "yes"){   
        #  $realORFstart=$exon_feature->start+($orf_start - 1);
        $realORFstart=$exon_feature->start+($orf_start - $total_exon_length_previous_round );
        my $end_part_of_exon=$exon_feature->start- $realORFstart + 1;
        if($end_part_of_exon >= $orf_length){           #exon      ============================================         
           $realORFend=$realORFstart+$orf_length-1;       #cds              =========================       
           last;
         }
        $mapped_length=$exon_feature->end - $realORFstart + 1;
        $first="no";
      }
      else{
        $mapped_length=$exon_feature->end - $exon_feature->start + 1;
      }
    }
    #exon are over the end of cds => we finish at this round
    if($total_exon_length >= ($orf_start+$orf_length) ){        #exon      ============================================  
      if($realORFstart > $exon_feature->start){                 #cds       ========================= 
        $realORFend=$realORFstart+$the_rest_to_map - 1 ;
      last;
      }else{
        $realORFend=$exon_feature->start + $the_rest_to_map - 1 ;
      last;
      }
    }
  }
return $realORFstart, $realORFend;
}

# Check the start and end of mRNA and gene feature;
sub check_start_end_of_mrna_feature{

  my ($mRNA_feature, $exon_list)=@_;

  ######
  #Modify mRNA start-end based on exon features
  my $exonStart=$exon_list->[0]->start;
  my $exonEnd=$exon_list->[$#{$exon_list}]->end;
  if ($mRNA_feature->start != $exonStart){
    $mRNA_feature->start($exonStart);
  }
  elsif($mRNA_feature->end != $exonEnd){
    $mRNA_feature->end($exonEnd);
  }
}

# Check the start and end of gene feature based on its mRNA;
sub check_start_end_of_gene_feature{

  my ($hash_omniscient, $gene_id)=@_;

  #####
  #Modify gene start-end (have to check size of each mRNA)
  my $geneExtremStart=1000000000000;
  my $geneExtremEnd=0;
  foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
    foreach my $mrna_feature ( @{$hash_omniscient->{'level2'}{'mrna'}{$gene_id}}) {
      my $start=$mrna_feature->start();
      my $end=$mrna_feature->end();

      if ($start < $geneExtremStart){
        $geneExtremStart=$start;
      }
      if($end > $geneExtremEnd){
        $geneExtremEnd=$end;
      }
    }
  }
  my $gene_feature=$hash_omniscient->{'level1'}{'gene'}{$gene_id};
  if ($gene_feature->start != $geneExtremStart){
      $gene_feature->start($geneExtremStart);
   }
   elsif($gene_feature->end != $geneExtremEnd){
      $gene_feature->end($geneExtremEnd);
    }
}


# The exons containing the original cds keep their parent names. The exon containing the new cds will have a new parent name.
sub create_two_exon_lists {
  # orignalFirst == true if original gene is first on the prediction
  my ($exons_features,$firstEnd, $secondStart, $orignalFirst)=@_;
  my @list_exon_originalPred;
  my @list_exon_newPred;

  foreach my $exon_feature (@$exons_features){ #for each exon
#    print "start:".$exon_feature->start." end:".$exon_feature->end."\n";
    if(two_positions_on_feature($exon_feature,$firstEnd,$secondStart)){ # We have to split the exon_feature P.S: We will loss sequence between the two positions
#      print "both on feature\n";
      my $duplicated_exon_feature=clone($exon_feature);#create a copy of the feature
      #manage original exon
      $exon_feature->end($firstEnd);
      $duplicated_exon_feature->start($secondStart);

      if($orignalFirst eq "true"){
        push( @list_exon_originalPred, $exon_feature);

        my @values = $duplicated_exon_feature->get_tag_values('ID');                 
        my $value = $values[0];
        create_or_replace_tag($duplicated_exon_feature,'ID', 'new_'.$value);  
        my @values = $duplicated_exon_feature->get_tag_values('Parent');                 
        my $value = $values[0];
        create_or_replace_tag($duplicated_exon_feature,'Parent', 'new_'.$value); 
        push( @list_exon_newPred, $duplicated_exon_feature);
        next;
      }else{ #original pred after
        push( @list_exon_originalPred, $duplicated_exon_feature);

        my @values = $exon_feature->get_tag_values('ID');                 
        my $value = $values[0];
        create_or_replace_tag($exon_feature,'ID', 'new_'.$value);  
        my @values = $exon_feature->get_tag_values('Parent');                 
        my $value = $values[0];
        create_or_replace_tag($exon_feature,'Parent', 'new_'.$value);
        push( @list_exon_newPred, $exon_feature);
        next;
      }
    }
    if(! (($exon_feature->end <=  $secondStart) and ($exon_feature->start >=  $firstEnd))){ # We remove it because exon between CDSs
      if ($exon_feature->end <=  $secondStart) { 
        if ($orignalFirst eq "true"){
          push( @list_exon_originalPred, $exon_feature);
        }else{
          my $duplicated_exon_feature=clone($exon_feature);#create a copy of the feature
          my @values = $duplicated_exon_feature->get_tag_values('ID');                 
          my $value = $values[0];
          create_or_replace_tag($duplicated_exon_feature,'ID', 'new_'.$value);
          my @values = $duplicated_exon_feature->get_tag_values('Parent');                 
          my $value = $values[0];
          create_or_replace_tag($duplicated_exon_feature,'Parent', 'new_'.$value);  
          push( @list_exon_newPred, $duplicated_exon_feature);
        }
      }
      if ($exon_feature->start >=  $firstEnd) { 
        if($orignalFirst eq "true"){
            my $duplicated_exon_feature=clone($exon_feature);#create a copy of the feature
            my @values = $duplicated_exon_feature->get_tag_values('ID');                 
            my $value = $values[0];
            create_or_replace_tag($duplicated_exon_feature,'ID', 'new_'.$value);  
            my @values = $duplicated_exon_feature->get_tag_values('Parent');                 
            my $value = $values[0];
            create_or_replace_tag($duplicated_exon_feature,'Parent', 'new_'.$value);  
            push( @list_exon_newPred, $duplicated_exon_feature);
          }else{
            push( @list_exon_originalPred, $exon_feature);
          }
        }
      }
  }
  my @list_exon_originalPred_sorted = sort {$a->start <=> $b->start} @list_exon_originalPred;
  my @list_exon_newPred_sorted = sort {$a->start <=> $b->start} @list_exon_newPred;

  return \@list_exon_originalPred_sorted, \@list_exon_newPred_sorted;
}

sub position_on_feature {

  my ($feature,$position)=@_;

  my $isOnSameExon=undef;
  if ( ($position >= $feature->start and $position <= $feature->end)){
    $isOnSameExon="true";
  }
  return $isOnSameExon;
}

sub two_positions_on_feature {

  my ($feature,$position1,$position2)=@_;

  my $areOnSameExon=undef;
  if ( ($position1 >= $feature->start and $position1 <= $feature->end) and ($position2 >= $feature->start and $position2 <= $feature->end) ){
    $areOnSameExon="true";
  }
  return $areOnSameExon;
}

sub translate_JD {
   my ($self,@args) = @_;
     my ($terminator, $unknown, $frame, $codonTableId, $complete,
     $complete_codons, $throw, $codonTable, $orf, $start_codon, $no_start_by_aa, $offset);

   ## new API with named parameters, post 1.5.1
   if ($args[0] && $args[0] =~ /^-[A-Z]+/i) {
         ($terminator, $unknown, $frame, $codonTableId, $complete,
         $complete_codons, $throw,$codonTable, $orf, $start_codon, $no_start_by_aa, $offset) =
       $self->_rearrange([qw(TERMINATOR
                                               UNKNOWN
                                               FRAME
                                               CODONTABLE_ID
                                               COMPLETE
                                               COMPLETE_CODONS
                                               THROW
                                               CODONTABLE
                                               ORF
                                               START
                                               NOSTARTBYAA
                                               OFFSET)], @args);
   ## old API, 1.5.1 and preceding versions
   } else {
     ($terminator, $unknown, $frame, $codonTableId,
      $complete, $throw, $codonTable, $offset) = @args;
   }
    
    ## Initialize termination codon, unknown codon, codon table id, frame
    $terminator = '*'    unless (defined($terminator) and $terminator ne '');
    $unknown = "X"       unless (defined($unknown) and $unknown ne '');
    $frame = 0           unless (defined($frame) and $frame ne '');
    $codonTableId = 1    unless (defined($codonTableId) and $codonTableId ne '');
    $complete_codons ||= $complete || 0;
    
    ## Get a CodonTable, error if custom CodonTable is invalid
    if ($codonTable) {
     $self->throw("Need a Bio::Tools::CodonTable object, not ". $codonTable)
      unless $codonTable->isa('Bio::Tools::CodonTable');
    } else {
        
        # shouldn't this be cached?  Seems wasteful to have a new instance
        # every time...
    $codonTable = Bio::Tools::CodonTable->new( -id => $codonTableId);
   }

    ## Error if alphabet is "protein"
    $self->throw("Can't translate an amino acid sequence.") if
    ($self->alphabet =~ /protein/i);

    ## Error if -start parameter isn't a valid codon
   if ($start_codon) {
     $self->throw("Invalid start codon: $start_codon.") if
      ( $start_codon !~ /^[A-Z]{3}$/i );
   }

   my $seq;

   if ($offset) {
    $self->throw("Offset must be 1, 2, or 3.") if
        ( $offset !~ /^[123]$/ );
    my ($start, $end) = ($offset, $self->length);
    ($seq) = $self->subseq($start, $end);
   } else {
    ($seq) = $self->seq();
   }

         ## ignore frame if an ORF is supposed to be found
   my $orf_region;
   if ( $orf ) {
            ($orf_region) = _find_orfs_nucleotide_JD( $self, $seq, $codonTable, $start_codon, $no_start_by_aa, $orf eq 'longest' ? 0 : 'first_only' );
            $seq = $self->_orf_sequence( $seq, $orf_region );
   } else {
   ## use frame, error if frame is not 0, 1 or 2
     $self->throw("Valid values for frame are 0, 1, or 2, not $frame.")
      unless ($frame == 0 or $frame == 1 or $frame == 2);
     $seq = substr($seq,$frame);
         }

    ## Translate it
    my $output = $codonTable->translate($seq, $complete_codons);
    # Use user-input terminator/unknown
    $output =~ s/\*/$terminator/g;
    $output =~ s/X/$unknown/g;

    ## Only if we are expecting to translate a complete coding region
    if ($complete) {
     my $id = $self->display_id;
     # remove the terminator character
     if( substr($output,-1,1) eq $terminator ) {
       chop $output;
     } else {
       $throw && $self->throw("Seq [$id]: Not using a valid terminator codon!");
       $self->warn("Seq [$id]: Not using a valid terminator codon!");
     }
     # test if there are terminator characters inside the protein sequence!
     if ($output =~ /\Q$terminator\E/) {
             $id ||= '';
       $throw && $self->throw("Seq [$id]: Terminator codon inside CDS!");
       $self->warn("Seq [$id]: Terminator codon inside CDS!");
     }
     # if the initiator codon is not ATG, the amino acid needs to be changed to M
     if ( substr($output,0,1) ne 'M' ) {
       if ($codonTable->is_start_codon(substr($seq, 0, 3)) ) {
         $output = 'M'. substr($output,1);
       }  elsif ($throw) {
         $self->throw("Seq [$id]: Not using a valid initiator codon!");
       } else {
         $self->warn("Seq [$id]: Not using a valid initiator codon!");
       }
     }
    }

    my $seqclass;
    if ($self->can_call_new()) {
     $seqclass = ref($self);
    } else {
     $seqclass = 'Bio::PrimarySeq';
     $self->_attempt_to_load_Seq();
    }
    my $out = $seqclass->new( '-seq' => $output,
                    '-display_id'  => $self->display_id,
                    '-accession_number' => $self->accession_number,
                    # is there anything wrong with retaining the
                    # description?
                    '-desc' => $self->desc(),
                    '-alphabet' => 'protein',
                              '-verbose' => $self->verbose
            );
    return $out, $orf_region;
}

sub concatenate_feature_list{

  my ($feature_list) = @_;

  my $seq = "";
  my $ExtremStart=1000000000000;
  my $ExtremEnd=0;

  foreach my $feature (@$feature_list) { 
#        my @values = $feature->get_tag_values('Parent');                 
#        my $parent = $values[0];
#        my @values = $feature->get_tag_values('Parent');                 
#        my $id = $values[0];
#    print $feature->primary_tag." ".$parent." ".$id."\n";
    my $start=$feature->start();
    my $end=$feature->end();
    my $seqid=$feature->seq_id();   
    $seq .= $db->seq( $seqid, $start, $end );

    if ($start < $ExtremStart){
      $ExtremStart=$start;
    }
    if($end > $ExtremEnd){
              $ExtremEnd=$end;
    }
  }
   return $ExtremStart, $seq, $ExtremEnd;
}

sub _find_orfs_nucleotide_JD {
    my ( $self, $sequence, $codon_table, $start_codon, $no_start_by_aa, $first_only ) = @_;
    $sequence    = uc $sequence;
    $start_codon = uc $start_codon if $start_codon;

    my $is_start = $start_codon
        ? sub { shift eq $start_codon }
        : sub { $codon_table->is_start_codon( shift ) };

    # stores the begin index of the currently-running ORF in each
    # reading frame
    my @current_orf_start = (-1,-1,-1);

    #< stores coordinates of longest observed orf (so far) in each
    #  reading frame
    my @orfs;

    # go through each base of the sequence, and each reading frame for each base
    my $seqlen = CORE::length $sequence;
    for( my $j = 0; $j <= $seqlen-3; $j++ ) {
        my $frame = $j % 3;

        my $this_codon = substr( $sequence, $j, 3 );
        my $AA = $codon_table->translate($this_codon);

        # if in an orf and this is either a stop codon or the last in-frame codon in the string
        if ( $current_orf_start[$frame] >= 0 ) {
            if ( _is_ter_codon_JD( $this_codon ) ||( my $is_last_codon_in_frame = ($j >= $seqlen-5)) ) {
                # record ORF start, end (half-open), length, and frame
                my @this_orf = ( $current_orf_start[$frame], $j+3, undef, $frame );
                my $this_orf_length = $this_orf[2] = ( $this_orf[1] - $this_orf[0] );

                $self->warn( "Translating partial ORF "
                                 .$self->_truncate_seq( $self->_orf_sequence( $sequence,\@ this_orf ))
                                 .' from end of nucleotide sequence'
                            )
                    if $first_only && $is_last_codon_in_frame;

                return\@ this_orf if $first_only;
                push @orfs,\@ this_orf;
                $current_orf_start[$frame] = -1;
            }
        }
        # if this is a start codon
        elsif ($is_start->($this_codon)) {
          if($no_start_by_aa){

            if($AA ne $no_start_by_aa){
              $current_orf_start[$frame] = $j;
            }
          }
          else{
            $current_orf_start[$frame] = $j;
          }
        }
    }

    return sort { $b->[2] <=> $a->[2] } @orfs;
}

# We can be sure it's a stop codon even with IUPAC
sub _is_ter_codon_JD{
  my ($codon) = @_;
  $codon=lc($codon);
  $codon =~ tr/u/t/;
  my $is_ter_codon=undef;

  if( ($codon eq 'tga') or ($codon eq 'taa') or ($codon eq 'tag') or ($codon eq 'tar') or ($codon eq 'tra') ){
    $is_ter_codon="yes";
  }
  return $is_ter_codon;
}

__END__

=head1 NAME

gff3_fixLongestORF.pl -
The script take a gff3 file as input. -
The script looks for other ORF in each gene model described in the gff file.
Several ouput files will be written if you specify an output. One will contain the gene not modified (intact), 
one the gene models fixed, one contains the putative pseudogene detected (As they are just putatuve, they are also present among the intacts ), and a last a report of the results.
Pseudogene particularity: If gene contains mRNA models goods and mRNA that look like a pseudogene, the pseudogene one will be removed. 


=head1 SYNOPSIS

    ./gff3_fixLongestORF.pl -gff=infile.gff --fasta genome.fa [ -o outfile ]
    ./gff3_fixLongestORF.pl --help

=head1 OPTIONS

=over 8

=item B<-gff>

Input GFF3 file that will be read (and sorted)

=item B<-fa> or B<--fasta>

Genome fasta file
The name of the fasta file containing the genome to work with.

=item B<-m> or B<--model>

Kind of ORF fix Model you want. By default all are used. To define specific model writte: --model 1,4
Model1= sequence original is part of new prediction; the predicted one is longest
Model2= sequence original predicted are different; the predicted one is longest, they don't overlap each other.
Model3= original protein and predicted one are different; the predicted one is longest, they overlap each other.
Model4= The prediction is shorter... /!\

=item B<-s> or B<--split>

This option is usefull for Model1 and Model2. Indeed when the prediction is non overlaping the original cds, it is possible to split the gene into two different genes. By default we don't split it.
We keep the longest. If you want to split it type: -s 

=item B<-o> , B<--output> , B<--out> or B<--outfile>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.

=back

=cut