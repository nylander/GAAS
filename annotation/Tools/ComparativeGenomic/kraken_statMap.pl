#!/usr/bin/env perl

use Carp;
use strict;
use File::Basename;
use Getopt::Long;
use Data::Dumper;
use Statistics::R;
use IO::File;
use Pod::Usage;
use List::MoreUtils qw(uniq);
use Clone 'clone';
use BILS::Handler::GFF3handler qw(:Ok);
use BILS::Handler::GXFhandler qw(:Ok);
use BILS::Plot::R qw(:Ok);

#####
# What we call parial gene (containing "_partial_part-" in the ID) ?
# This gene has been seen as patial: During a lift-over a gene can be detected on 2 several contigs.
#                                    (full kraken file => features with kraken attribute to TRUE are on contig of the target genome (Transfert annotation on), the others (kraken attribute to FALSE) are on the reference genome to liftfover (where annotations are taken to try to liftover) )
#####


my $usage = qq{
########################################################
# BILS 2015 - Sweden                                   #  
# jacques.dainat\@bils.se                               #
# Please cite BILS (www.bils.se) when using this tool. #
########################################################
};

my $outfile = undef;
my $gff = undef;
my $valueK = undef;
my $attributes = undef ;
my $help= 0;

if ( !GetOptions(
    "help|h" => \$help,
    "gff=s" => \$gff,
    "value|v|threshold|t=i" => \$valueK,
    "outfile|output|out|o=s" => \$outfile))

{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($help) {
    pod2usage( { -verbose => 2,
                 -exitval => 0,
                 -message => "$usage\n" } );
}
 
if ( ! (defined($gff)) ){
    pod2usage( {
           -message => "\nAt least 1 parameter is mandatory:\nInput reference gtf file (--f)\n\n".
           "$usage\n",
           -verbose => 0,
           -exitval => 2 } );
}


## Manage output file
$outfile=~ s/.gff//g;

my $ostreamPlotFile = new IO::File;
my $pathPlotFile=$outfile."-geneMapped.txt";
my $pathOutPlot=$outfile."-geneMapped_plot.pdf"; 
$ostreamPlotFile->open($pathPlotFile, 'w' ) or
        croak(
          sprintf( "Can not open '%s' for writing %s", $pathPlotFile, $! )
        );

## Manage output file
my $gffout;
my $outReport     = IO::File->new();
if ($outfile) {
  open(my $fh, '>', $outfile.".gff") or die "Could not open file '$outfile' $!";
  $gffout= Bio::Tools::GFF->new(-fh => $fh, -gff_version => 3 );

  $outReport->open($outfile."_report.txt", 'w') or die "Could not open file '$outfile'_report.txt $!";
}
else{
  $gffout = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 3);

  $outReport->fdopen( fileno(STDOUT), 'w' ) or die "Could not open file STDOUT $!";
}

my $messageValue;
if ( defined($valueK) ){
  $messageValue = "You choose to keep in output only genes mapped over $valueK percent.\n" 
}else{
  $messageValue = "We will keep all the mapped features.\n";
  $valueK=0;
}

#print info
if ($outfile) {
  print $outReport $messageValue;
  }else{print $messageValue;}

                #####################
                #     MAIN          #
                #####################


######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) = BILS::Handler::GXFhandler->slurp_gff3_file_JD($gff);


my $nbOriginalGene = nb_feature_level1($hash_omniscient);
my %mappedPercentPerGene;
my %seqID;
my @new_list_l3;
my $nb_split=0;
my $kraken_tag = undef;

#####
# == LEVEL 1
######
foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
  foreach my $id_tag_key_level1 (keys %{$hash_omniscient->{'level1'}{$primary_tag_key_level1}}){
    
    my $gene_feature = $hash_omniscient->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1};
    my @ListmrnaNoMatch;
    my $at_least_one_level2_match="no";
    
    #####
    # == LEVEL 2
    ######
    my $new_l2_seqID=undef;
    foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
        
      if ( exists_keys($hash_omniscient, ('level2',$primary_tag_key_level2,$id_tag_key_level1) ) ){
        foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}) {
          
          my $percentMatch=0;
          my $totalSize=0;
          my $matchSize=0;

          my $level2_ID = lc($feature_level2->_tag_value('ID'));

          ######
          # == LEVEL 3
          # We will look the size of mapped features
          # Feature can be exon ,cds or utr in that order 
          my $refListFetaureL3 = takeOneListLevel3From1idLevel2($hash_omniscient, $level2_ID);
          if(! $refListFetaureL3 ){next;}

          my $atLeastOneMapped=undef;
          
          undef %seqID;
          my %l3LengthBySeqID=undef;
          foreach my $feature (@{$refListFetaureL3}){  
            my $end=$feature->end();
            my $start=$feature->start();

            if(! $kraken_tag ){
              if($feature->has_tag('kraken_mapped')){
                  $kraken_tag = "kraken_mapped";
              }
              elsif($feature->has_tag('Kraken_mapped')){
                  $kraken_tag = "Kraken_mapped";
              }
            }

            my $mapping_state = undef;
            if($feature->has_tag($kraken_tag)){
              $mapping_state = lc($feature->_tag_value($kraken_tag));
            }
            else{ print "error !! No $kraken_tag attribute found for the feature".$feature->gff_string()."\n";}


            if( $mapping_state eq "true"){
                $seqID{$feature->seq_id()}++;
                $l3LengthBySeqID{$feature->seq_id()}+=($end-$start);
                $atLeastOneMapped=1;
            }
            elsif(! $mapping_state eq "false"){
              print "error !! We don't understand the $kraken_tag attribute value found for the feature".$feature->gff_string()."\n Indeed, we expect false or true.\n";
            }

            $totalSize+=($end-$start);   
          }

            
          if(! $atLeastOneMapped){ #we have to remove the mRNA
            push (@ListmrnaNoMatch, $level2_ID);
            next; # next mRNA
          }

          #in case there is several seqID
          foreach my $seqIDkey (keys (%l3LengthBySeqID)){

            my $matchSize = $l3LengthBySeqID{$seqIDkey};
            #compute the MATCH
            $percentMatch=($matchSize*100)/$totalSize;
            
            #######
            # Add information to gff
            ########
            if ($percentMatch >= $valueK) {   
              #We print gene only if a percentage match value is superior to the threshold fixed (No threshold equal everything = 0)
              $at_least_one_level2_match="yes";
              $percentMatch = sprintf('%.2f', $percentMatch);

              manage_gene_label($gene_feature, $percentMatch, $kraken_tag); # add info to level1 (gene) feature 

              create_or_replace_tag($feature_level2,$kraken_tag,$percentMatch."%"); # add info to level2 (mRNA) feature
              create_or_replace_tag($feature_level2,'description',"Mapped at ".$percentMatch."%"); # add info to level2 (mRNA) feature

              #save best value for gene
              if(!exists($mappedPercentPerGene{$id_tag_key_level1})){ # case where it doesn t exist
                $mappedPercentPerGene{$id_tag_key_level1}=$percentMatch;
              }
              elsif($mappedPercentPerGene{$id_tag_key_level1} < $percentMatch){ # case where it exists but better value to save
                $mappedPercentPerGene{$id_tag_key_level1}=$percentMatch;
              }
            }
          }

          #Now remove element form L3 list that have not been mapped
          my @id_concern_list=($level2_ID);
          my @tag_list=('all');
          remove_element_from_omniscient_attributeValueBased(\@id_concern_list, 'false', $kraken_tag, $hash_omniscient, 'level3', 'false', \@tag_list);

          #check SeqID
          if(! exists($seqID{$feature_level2->seq_id()})){ #If the seqid of the level1 was the one from the ref
            foreach my $key (sort(keys %seqID)){
              $feature_level2->seq_id($key);
              $new_l2_seqID=$key;
            
            }
          }

          my $nb_seqID = keys %seqID;
          if($nb_seqID > 1){
            $nb_split++;
            foreach my $seqID (keys %seqID){
              if($seqID ne $feature_level2->seq_id()){

                ###########
                #
                foreach my $tag_l3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
                    my @ListIDtoRemove=();
                    if ( exists ($hash_omniscient->{'level3'}{$tag_l3}{$level2_ID} ) ){
                      foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{$tag_l3}{$level2_ID}}) {
                        if($feature_level3->seq_id eq $seqID){

                          my $new_parent=$feature_level3->seq_id().$feature_level3->_tag_value('Parent'); # we need a new parent to deal the reconstrunction of correct features
                          create_or_replace_tag($feature_level3, 'Parent', $new_parent);
                          push @new_list_l3, clone($feature_level3);
                          push @ListIDtoRemove, $feature_level3->_tag_value('ID');
                        }
                      }
                    }
                    # Remove from this hash element that are not mapped with the main stream
                    my @ListGeneID=($level2_ID);
                    my @listTag=($tag_l3);
                    remove_element_from_omniscient(\@ListGeneID, \@ListIDtoRemove, $hash_omniscient, 'level3', 'true', \@listTag);
                }
              }
            }
          }
        }
      }
    }
    ##########
    # all mRNA of the gene has been tested and we don't have any match over the choosen treshold
    if ($at_least_one_level2_match eq "no"){ 
      delete $hash_omniscient->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1}; # remove gene under threshold fixed
    }
    ##################
    #at least one mRNA with at least one exon (or cds or utr) mapped. 
    else{
      #check the gene is the same seqID as L2 if L2 has been modified
      if($new_l2_seqID and $gene_feature->seq_id ne $new_l2_seqID){
        $gene_feature->seq_id($new_l2_seqID);
      }

      #remove non-mapped L2
      if($#ListmrnaNoMatch != -1 ){ #one mRNA didnt match we have to remove them
              my @ListGeneID=($id_tag_key_level1);
              my @listTag=('mrna');
              remove_element_from_omniscient(\@ListGeneID, \@ListmrnaNoMatch, $hash_omniscient, 'level2', 'true', \@listTag);
      }
    }
  }
}
#Add to the hash amn extra list pf feature to take care of
$hash_omniscient->{'extra'} = \@new_list_l3;
# reinitialyzer the hash
($hash_omniscient, $hash_mRNAGeneLink) = BILS::Handler::GXFhandler->slurp_gff3_file_JD($hash_omniscient);
print "Calcul of mapped percentage length finished !\n";


######################
# Check if nothing mapped
my $nbKey = keys %mappedPercentPerGene;
if ($nbKey == 0){
 print "No succefully mapped feature found!\n"; exit;
}

###############
# print the value per gene in a temporary file for R plot
foreach my $key (keys %mappedPercentPerGene){
   print $ostreamPlotFile "$mappedPercentPerGene{$key}\n";
}


########
#print GFF the selected features (over the choosen treshold)
########
print_omniscient($hash_omniscient, $gffout);

my $nbGeneMapped= keys %mappedPercentPerGene;

my $messageEnd;
$messageEnd.= "\nTo resume:\n==========\n\n";
$messageEnd.= "The original file contained $nbOriginalGene genes\n";
$messageEnd.= "$nbGeneMapped of them have been mapped over the $valueK % match threshold. \n";
my $nbEndGene = nb_feature_level1($hash_omniscient);
if($nbEndGene != $nbOriginalGene){
  my $nb_split = $nbEndGene -$nbOriginalGene;
  $messageEnd.= "Among them, we have $nb_split cases where the gene was mapped at several places (at least two different scaffold/Contig/chromosome). In total we have $nb_split split \n";
}
$messageEnd.= "\n";

#print info
if ($outfile) {
print $outReport $messageEnd;
}else{print $messageEnd;}

#############
#PLOT
#############
# Create the legend
my $nbOfGeneSelected= $nbGeneMapped;
# parse file name to remove extension
my ($file1,$dir1,$ext1) = fileparse($gff, qr/\.[^.]*/);
my $legend=$nbOfGeneSelected." genes selected from ".$file1;

my @listTuple=([$pathPlotFile,$legend]);
my $R_command=rcc_density_one_row_per_file(\@listTuple,"histogram","Percentage of gene length mapped","10","",$pathOutPlot); # create the R command
execute_R_command($R_command);

my $messagePlot;
$messagePlot = "Plot done in the pdf file named $pathOutPlot\n";

#print info
if ($outfile) {
  print $outReport $messagePlot;
}else{print $messagePlot;}

# Delete temporary file
unlink "$pathPlotFile";

#END
print "We finished !! Bye Bye.\n";


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


# Feature we look at are in the order exon ,cds or utr 
#
sub takeOneListLevel3From1idLevel2 {
  my ($hash_omniscient, $level2_ID)=@_;
  
  my  $refListFetaureL3=undef;
 
  if ( exists_keys($hash_omniscient, ('level3','exon',$level2_ID) ) ){ 
    $refListFetaureL3 = $hash_omniscient->{'level3'}{'exon'}{$level2_ID};
  }
  elsif ( exists_keys($hash_omniscient,('level3','cds',$level2_ID) ) ){ 
    $refListFetaureL3 = $hash_omniscient->{'level3'}{'cds'}{$level2_ID};
  }
  else{
    my $match=undef;
    foreach my $tag (keys %{$hash_omniscient->{'level3'}}){
      if($hash_omniscient->{'level3'}{$tag}{$level2_ID}){
        $match="yes";
        if($tag =~ "utr"){
          $refListFetaureL3 = $hash_omniscient->{'level3'}{$tag}{$level2_ID};
        }
        #else{}
      }
    }
    if(! $match){
      print "No feature level3 expected found for ".$level2_ID." level2 !\n";
    }
  }
  return  $refListFetaureL3;
}

sub manage_gene_label{

  my ($gene_feature, $percentMatch, $kraken_tag)=@_;
  if (! $gene_feature->has_tag($kraken_tag)){ # No kraken_mapped attribute
    label_by_value($gene_feature, $percentMatch, $kraken_tag);
  }
  else{ # kraken_mapped tag exists, check if we have to change it
    my @values = $gene_feature->get_tag_values($kraken_tag);
    my $alreadyMap = lc(shift @values) ;
    if ($alreadyMap eq "false" or $alreadyMap eq "true"){ 
      label_by_value($gene_feature, $percentMatch, $kraken_tag);
    }
    elsif ( ($alreadyMap ne 'full') and (( $percentMatch != 0 ) and ($alreadyMap eq 'none')) ){ # if the existing tag is full or the new tag we want to add is none ($percentMatch == 0), we skip it.
      create_or_replace_tag($gene_feature,$kraken_tag,'partial'); # add info to gene feature
    }
  }
}

sub label_by_value{
  my ($gene_feature, $percentMatch, $kraken_tag)=@_;

  if($percentMatch == 100){
    create_or_replace_tag($gene_feature,$kraken_tag,'full'); # add info to gene feature
  }
  elsif ($percentMatch != 0){
    create_or_replace_tag($gene_feature,$kraken_tag,'partial'); # add info to gene feature
  }
  else{
    create_or_replace_tag($gene_feature,$kraken_tag,'none'); # add info to gene feature
  }
}

sub plotR {
  my ($pathIn, $pathOut, $title)=@_;

  my $R = Statistics::R->new() or die "Problem with R : $!\n";

#R command
$R->send(
      qq`
      listValues=as.matrix(read.table("$pathIn", sep="\t", he=F))
      legendToDisplay=paste("Number of value used : ",length(listValues))

      pdf("$pathOut")
      hist(listValues, breaks=seq(0,100,5), xlab="Percentage of cds lift-over", ylab="Number", main="$title")
      dev.off()`
          );

# Close the bridge
$R->stopR();

# Delete temporary file
unlink "$pathIn";
}


__END__

=head1 NAME

kraken_statMap.pl -
The script take a gff file as input. It will analyse the kraken_mapped attributes to calculate the mapped percentage of each mRNA. 
/!\ The script handles chimeric files (i.e containg gene part mapped on the template genome and others on the de-novo one) 
/!\/¡\ If the file is complete (containing kraken_mapped="TRUE" and kraken_mapped="FALSE" attributes), the script calcul the real percentage lentgh that has been mapped. Else the calcul is only based on feature with kraken_mapped="TRUE" attributes. So in this case the result most of time will be 100%, execpt for cases where several piecies are mapped at different area of the de-novo genome.
According to a threshold (0 by default), gene with a mapping percentage over that value will be reported.
A plot nammed geneMapped_plot.pdf is performed to visualize the result. 

=head1 SYNOPSIS

    ./kraken_statMap.pl -gtf=infile.gff [ -o outfile ]
    ./kraken_statMap.pl --help

=head1 OPTIONS

=over 8

=item B<-gff>

Input gff file that will be read.

=item B<-v>, B<--value>, B<--threshold> or B<-t>

Gene mapping percentage over which a gene must be reported. By default the value is 0.

=item B<-o> , B<--output> , B<--out> or B<--outfile>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.

=back

=cut
