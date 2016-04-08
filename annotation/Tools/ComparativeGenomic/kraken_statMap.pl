#!/usr/bin/perl

use Carp;
use strict;
use File::Basename;
use Getopt::Long;
use Data::Dumper;
use Statistics::R;
use IO::File;
use Pod::Usage;
use List::MoreUtils qw(uniq);
use Bio::Tools::GFF;
use BILS::Handler::GXFhandler qw(:Ok);
use BILS::Handler::GFF3handler qw(:Ok);
use BILS::Plot::R qw(:Ok);

#####
# What we call parial gene (containing "_partial_part-" in the ID) ?
# This gene has been seen as patial: During a lift-over a gene can be detected on 2 several contigs. In the two different part are correctly labelled as partial
#                                    If the file in input is a chimera (full kraken file => features with kraken attribute to TRUE are on contig of the reference genome (Transfert annotation on), the others (kraken attribute to FALSE) are on the genome to liftfover (where annotations are taken to try to liftover) )
#                                       In that case this is a false partial gene. Because the gebe has not been mapped on several contig of the genome.
# In this script we correct these name details.
#####


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
		The name of the gff3 file to work with. >>> It must contain kraken_mapped attribute. <<<
		
  Ouput:    
    [--out filename]
        The name of the output file (A GFF file).

  Option:
    [--value Integer]
        Only mapped gene over that value will be reported on the output. (Exon features are considered first, If none then we try to work with CDS and still none we look for UTR.).

  /!\\ If Kraken_mapped attributes exist for gene or mRNA, they will be overwritte. The value between gene and mRNA may change (gene will keep only the higher isoform mapped value) !
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
my ($hash_omniscient, $hash_mRNAGeneLink) = BILS::Handler::GFF3handler->slurp_gff3_file_JD($gff);

##########
# merge same original mRNA together (was part of the same gene before kraken)
# Needed to know the total size of the mRNA before to have been split (name modified) because don't map or map on other Contig during kraken process.
#########
my %OldCurrentGeneName;
my %CheckIfReallyPartialMap;
foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
  foreach my $id_tag_key_level1 (keys %{$hash_omniscient->{'level1'}{$primary_tag_key_level1}}){
     
    foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
      foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}) {

        my @values = $feature_level2->get_tag_values('ID');
        my $current_mRNA_name = lc(shift @values) ;

        my @originalName_list = split /_partial_part-/, $current_mRNA_name;
        my $originalName=$originalName_list[0];
        push(@{$OldCurrentGeneName{$originalName}},$current_mRNA_name);
      }
    }
  }
}



my %mappedPercentPerGene;
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
    foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
        
      if ( exists ($hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1} ) ){
        foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}) {
          
          my $percentMatch=0;
          my $totalSize=0;
          my $matchSize=0;

          my @temp = $feature_level2->get_tag_values('ID');
          my $level2_ID_original= shift @temp;
          my $level2_ID = lc($level2_ID_original);

          ######
          # == LEVEL 3
          # We will look the size of mapped features
          # Feature can be exon ,cds or utr in that order 
          my $refListFetaureL3=takeOneListLevel3From1idLevel2($hash_omniscient, $level2_ID);
          if(! $refListFetaureL3 ){next;}

          my $atLeastOneMapped=undef;
          foreach my $exon_feature (@$refListFetaureL3){  
            my $end=$exon_feature->end();
            my $start=$exon_feature->start();

            my @values = $exon_feature->get_tag_values('kraken_mapped');
            if($#values == -1){print "error !! No kraken_mapped attribute found for the feature $exon_feature->gff_string()\n";}
            my $mapping_state = lc(shift @values) ;

            if( $mapping_state eq "true"){
                $matchSize+=($end-$start);
                $atLeastOneMapped=1;
            }
            $totalSize+=($end-$start);   
          }
            
          if(! $atLeastOneMapped){ #we have to remove the mRNA
            push (@ListmrnaNoMatch, $level2_ID_original);
            next; # next mRNA
          }

          ##########
          # Calcul the original size of the mRNA if have been split during the kraken process
          #########
          my @originalName_list = split /_partial_part-/, $level2_ID;
          my $originalName=$originalName_list[0];
          my $nbPieces = $#{$OldCurrentGeneName{$originalName}};
          if ( $nbPieces > 0 ){ # counter part gene exist (Maybe still on original genome)
            foreach my $mRNA_counterpart (@{$OldCurrentGeneName{$originalName}}){
              if ($mRNA_counterpart ne $level2_ID){
                #Take the correct list of feature
                my $refListFetaureL3=takeOneListLevel3From1idLevel2($hash_omniscient, $mRNA_counterpart);
                if(! $refListFetaureL3 ){next;}
                foreach my $exon_feature (@$refListFetaureL3){# exon ,cds or utr in that order        
                  my $end=$exon_feature->end();
                  my $start=$exon_feature->start();
                  $totalSize+=($end-$start);   
                }
              }
            }
          }
          $percentMatch=($matchSize*100)/$totalSize;
          
          #######
          # Add information to gff
          ########
          if ($percentMatch >= $valueK) {   
            #We print gene only if a percentage match value is superior to the threshold fixed (No threshold equal everything = 0)
            $at_least_one_level2_match="yes";
            $percentMatch = sprintf('%.2f', $percentMatch);

            manage_gene_label($gene_feature,$percentMatch); # add info to level1 (gene) feature 

            create_or_replace_tag($feature_level2,'kraken_mapped',$percentMatch."%"); # add info to level2 (mRNA) feature
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
      # We save the short gene name (originalName)!
      # Look if mapped gene is really _partial_part- (should be _partial_part- only if mapped TRUE on 2 different area of the genome). So if value=2 or more it's a real partial one.
      my @values = $gene_feature->get_tag_values('ID');
      my $Full_gene_name = lc(shift @values) ;

      my @originalName_list = split /_partial_part-/, $Full_gene_name;
      my $originalName=$originalName_list[0];
      $CheckIfReallyPartialMap{$originalName}++;

      if($#ListmrnaNoMatch != -1 ){ #one mRNA didnt match we have to remove them
              my @ListGeneID=($Full_gene_name);
              my @listTag=('mrna');
              remove_element_from_omniscient(\@ListGeneID, \@ListmrnaNoMatch, $hash_omniscient, 'level2', 'true', \@listTag);
      }
    }
  }
}
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


######################################################
# Manage ID and Parent attributes for non real partial mapped gene
my %info_split;
foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
  foreach my $id_tag_key_level1 (keys %{$hash_omniscient->{'level1'}{$primary_tag_key_level1}}){

    my $gene_feature = $hash_omniscient->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1};
    my @values = $gene_feature->get_tag_values('ID');
    my $current_gene_id = shift @values ;
    my @originalName_list = split /_partial_part-/, $current_gene_id;
    my $originalID=$originalName_list[0];

    #save some information
    if(($CheckIfReallyPartialMap{lc($originalID)} >= 2) and ($current_gene_id =~ m/_partial_part-/)){
      $info_split{$CheckIfReallyPartialMap{lc($originalID)}}++;
    }

    # We have to remove the partial-part information because it's wrong. It's due it was on two different genomes and not two different area of the same genome.
    if(($CheckIfReallyPartialMap{lc($originalID)} < 2) and ($current_gene_id =~ m/_partial_part-/)){

      create_or_replace_tag($gene_feature,'ID',$originalID); 
    
      foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
        foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}) {
          create_or_replace_tag($feature_level2,'Parent',$originalID);

          my @values_level2 = $feature_level2->get_tag_values('ID');
          my $current_mrna_id = shift @values_level2 ;
          my @originalmRNA_list = split /_partial_part-/, $current_mrna_id;
          my $original_mRNA_ID=$originalmRNA_list[0];
          create_or_replace_tag($feature_level2,'ID',$original_mRNA_ID);

          foreach my $primary_tag_key_level3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
            foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{$primary_tag_key_level3}{lc($current_mrna_id)}}) {
              create_or_replace_tag($feature_level3,'Parent',$original_mRNA_ID);

              my @values_level3 = $feature_level3->get_tag_values('ID');
              my $exon_id = shift @values_level3 ;
              $exon_id=~ s/_partial_part-//g;
              create_or_replace_tag($feature_level3,'ID',$exon_id);

            }
            @{$hash_omniscient->{'level3'}{$primary_tag_key_level3}{lc($original_mRNA_ID)}} = @{$hash_omniscient->{'level3'}{$primary_tag_key_level3}{lc($current_mrna_id)}}; # As to print level3 we based our search on ID of level2 we have to change that
            delete  $hash_omniscient->{'level3'}{$primary_tag_key_level3}{lc($current_mrna_id)};# delete ancient key level3
          }
        }
        @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$originalID}} = @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}; # As to print level2 we based our search on ID of level1 we have to change that
        delete  $hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1};# delete ancient key level2
      }
    $hash_omniscient->{'level1'}{$primary_tag_key_level1}{$originalID} = $hash_omniscient->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1}; # As to print level2 we based our search on ID of level1 we have to change that
    delete $hash_omniscient->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1};# delete ancient key level1
    }
  }
}
print "Manage ID and Parent attributes for non real partial mapped gene finished!\n";

######################################
# Check seq_id, start and stop for each feature level1 and level2 (Case of a full kraken file done with a gtf3)
# Filtering of level3 features. Must contain kraken_mapped="true";

foreach my $tag_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
  foreach my $id_level1 (keys %{$hash_omniscient->{'level1'}{$tag_level1}}){
    my $gene_feature=$hash_omniscient->{'level1'}{$tag_level1}{$id_level1};
    my $NeedFixIt=undef;           
    foreach my $tag_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level1 = gene or repeat etc...
      if ( exists ($hash_omniscient->{'level2'}{$tag_level2}{$id_level1} ) ){
        foreach my $feature_level2 (@{$hash_omniscient->{'level2'}{$tag_level2}{$id_level1}}){
          
          my @values = $feature_level2->get_tag_values('ID');
          my $level2_ID = lc(shift @values) ;

          foreach my $tag_level3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
            if ( exists ($hash_omniscient->{'level3'}{$tag_level3}{$level2_ID} ) ){
              my @NewListLevel3;
              foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{$tag_level3}{$level2_ID}}) {

                my @values = $feature_level3->get_tag_values('kraken_mapped');
                my $mapping_state = lc(shift @values) ;
                if( $mapping_state eq "true"){
                  if(($feature_level3->seq_id ne $gene_feature->seq_id) and !($NeedFixIt)){
                    $NeedFixIt=$feature_level3->seq_id;
                  }
                  push(@NewListLevel3 ,$feature_level3);
                }
              }
              if ($NeedFixIt){
                @{$hash_omniscient->{'level3'}{$tag_level3}{$level2_ID}} = @NewListLevel3;
              }
            }
          }
        if(($feature_level2->seq_id ne $NeedFixIt) and ($NeedFixIt)){
          $feature_level2->seq_id($NeedFixIt);
         
          # check start and stop
          check_level2_positions($hash_omniscient, $feature_level2);
        }
        }
      }
    }
  if($NeedFixIt){
    $gene_feature->seq_id($NeedFixIt);
   
    # check start and stop
    check_level1_positions($hash_omniscient,$gene_feature)
  }
  }
}
print "The checking of seq_id, start and stop information of each feature level1 and level2 finished!\n";

########
#print GFF the selected features (over the choosen treshold)
########
print_omniscient($hash_omniscient, $gffout);

my $nbGeneMapped= keys %mappedPercentPerGene;
my $nbOriginalGene=keys %OldCurrentGeneName;

my $messageEnd;
$messageEnd.= "\nTo resume:\n==========\n\n";
$messageEnd.= "The original file contained $nbOriginalGene genes\n";
$messageEnd.= "$nbGeneMapped of them have been mapped over the $valueK % match threshold. \n";
foreach my $nb_split (keys %info_split){
  $messageEnd.= "Among them, $info_split{$nb_split} genes have been mapped $nb_split times (2 different scaffold/Contig/chromosome)\n";
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

sub takeOneListLevel3From1idLevel2 {
  my ($hash_omniscient, $level2_ID)=@_;
  
  my  $refListFetaureL3=undef;
 
  if ( exists ($hash_omniscient->{'level3'}{'exon'}{$level2_ID} ) ){ 
    $refListFetaureL3 = $hash_omniscient->{'level3'}{'exon'}{$level2_ID};
  }
  elsif ( exists ($hash_omniscient->{'level3'}{'cds'}{$level2_ID} ) ){ 
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

  my ($gene_feature, $percentMatch)=@_;
  if (! $gene_feature->has_tag('kraken_mapped')){ # No kraken_mapped attribute
    label_by_value($gene_feature, $percentMatch);
  }
  else{ # kraken_mapped tag exists, check if we have to change it
    my @values = $gene_feature->get_tag_values('kraken_mapped');
    my $alreadyMap = lc(shift @values) ;
    if ($alreadyMap eq "false" or $alreadyMap eq "true"){ 
      label_by_value($gene_feature, $percentMatch);
    }
    elsif ( ($alreadyMap ne 'full') and (( $percentMatch != 0 ) and ($alreadyMap eq 'none')) ){ # if the existing tag is full or the new tag we want to add is none ($percentMatch == 0), we skip it.
      create_or_replace_tag($gene_feature,'kraken_mapped','partial'); # add info to gene feature
    }
  }
}

sub label_by_value{

  my ($gene_feature, $percentMatch)=@_;
  if($percentMatch == 100){
    create_or_replace_tag($gene_feature,'kraken_mapped','full'); # add info to gene feature
  }
  elsif ($percentMatch != 0){
    create_or_replace_tag($gene_feature,'kraken_mapped','partial'); # add info to gene feature
  }
  else{
    create_or_replace_tag($gene_feature,'kraken_mapped','none'); # add info to gene feature
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
Accodring to a threshold (0 by default), gene with a mapping percentage over that value will be reported.
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