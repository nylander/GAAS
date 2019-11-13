#!/usr/bin/env perl

## TO DO:
## Need to build UTR features from the difference
## between CDS and exon features.
use Carp;
use strict;
use Getopt::Long;
use Pod::Usage;
use List::MoreUtils qw(uniq);
use lib $ENV{ANDREASCODE};
use Private::Bio::IO::GFF;

my $usage = qq{
########################################################
# NBIS 2015 - Sweden                                   #	
# jacques.dainat\@nbis.se                               #
# Please cite NBIS (www.nbis.se) when using this tool. #
########################################################

That script lift-over the transcript GO term to the gene level. (duplicates are removed)

Usage: perl my_script.pl --gff Infile [--out outfile]
  Getting help:
    [--help]

  Input:
    [--gff filename]
		The name of the gff3 file to work with. 
		
  Ouput:    
    [--out filename]
        The name of the output file (A GFF file).

};

my $outfile = undef;
my $gff = undef;
my $valueK = undef;
my $attributes = undef ;
my $help= 0;

if ( !GetOptions(
    "help|h" => \$help,
    "gff=s" => \$gff,
    "value|v=i" => \$valueK,
    "outfile|out|o=s" => \$outfile))

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

if ( defined($valueK) ){
  print "You choose to keep in output only genes mapped over $valueK percent.\n" 
}

## Manage output file
my $ostreamPlotFile = IO::File->new();

my $pathPlotFile="geneMapped.txt";
my $pathOutPlot="geneMapped_plot.pdf"; 
$ostreamPlotFile->open($pathPlotFile, 'w' ) or
        croak(
          sprintf( "Can not open '%s' for writing %s", $pathPlotFile, $! )
        );

my $ostream = IO::File->new();
if ($outfile) {
  $ostream->open( $outfile, 'w' ) or
    croak( sprintf( "Can not open '%s' for writing %s", $outfile, $! ) );
}
else{
  $ostream->fdopen( fileno(STDOUT), 'w' ) or
    croak( sprintf( "Can not open STDOUT for writing: %s", $! ) );
}

my $output = Private::Bio::IO::GFF->new( ostream => $ostream );

### Parse GFF input file and add annotations
# Manage input gff3 file
my $ref_istream = IO::File->new();
$ref_istream->open( $gff, 'r' ) or
  croak(
     sprintf( "Can not open '%s' for reading: %s", $gff, $! ) );
my $ref_in = Private::Bio::IO::GFF->new(istream => $ref_istream);
# declaration hashes reference
my $ref_genes; my $refmRNA;my $refexon; my $refcds; my $refUTR5; my $refUTR3; my $refUTR; my $reftRNA; my $refRepeat; my $refRepeatMatch_part; my $refpieceStudied;


#Parse GFF 
($ref_genes,$refmRNA,$refexon, $refcds, $refUTR5, $refUTR3, $refUTR, $reftRNA, $refRepeat, $refRepeatMatch_part, $refpieceStudied) = parseGFF ($ref_in, $gff);
print "Parsing FINISH\n";

foreach my $geneName (keys %$ref_genes){ # For each gene
 printGene($geneName);
}
  

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


# method dedicated to the Exon,UTR,CDS printing
sub brickToprintHash{
  my ($featureList)=@_;
  foreach my $feature (@$featureList){
        #print feature
        $output->write_feature($feature);
  }
}

# method allowing to print all gene information (mRNAs,CDSs,Exons,UTRs...)
sub printGene {
  my ($geneName) = @_;

  # print gene 
  my $geneFeature=$ref_genes->{$geneName};
    #get GO term
    my @listGO;
    my @transcriptList=@{$refmRNA->{$geneName}};
    foreach my $transcriptFeature (@transcriptList){ # For each transcript
      $geneFeature->add_attribute('Ontology_term',$transcriptFeature->get_attribute('Ontology_term'));
    }
    $output->write_feature($geneFeature);
  
  # print mRNA
  foreach my $transcriptFeature (@transcriptList){ # For each transcript 
    $output->write_feature($transcriptFeature);

    #print other feature (CDS,UTR,EXON)
    my $transcriptName=$transcriptFeature->get_attribute('ID');
      #Get exonhash of the current mRNA studied
    if(exists $refexon->{$transcriptName}){      # If exon are compacted some mRNA could not have exons ...
      my @exonList=@{$refexon->{$transcriptName}};   
      brickToprintHash(\@exonList);
    }
    #Get cdshash of the current mRNA studied
    if (exists $refcds->{$transcriptName}){
      my @cdsList=@{$refcds->{$transcriptName}};
      brickToprintHash(\@cdsList);
    }
    #Get UTR5hash of the current mRNA studied if exist
    if (exists $refUTR5->{$transcriptName}){
      my @UTR5List=@{$refUTR5->{$transcriptName}};
      brickToprintHash(\@UTR5List);
    }
    #Get UTR3hash of the current mRNA studied if exist
    if (exists $refUTR3->{$transcriptName}){
      my @UTR3List=@{$refUTR3->{$transcriptName}};
      brickToprintHash(\@UTR3List);
    }
    #Get UTR3hash of the current mRNA studied if exist
    if (exists $refUTR->{$transcriptName}){
      my @UTRHash=@{$refUTR->{$transcriptName}};
      brickToprintHash(\@UTRHash);
    }
  }
}

sub addDataToHashOfHash {
  my ($hashOfHash,$key,$data) = @_;

  #Test if exon hash already exists
  if (exists $hashOfHash->{$key}){
    # put hash in hash
    push (@{$hashOfHash->{$key}} ,  $data);
  }
  else{          
    $hashOfHash->{$key} =  [$data];
  }
}

# method to parse GFF3 files
# take in account features gens,mRNA,tRNA,exon,CDS,three_prime_UTR and five_prime_UTR
sub parseGFF {
  my($file_in,$fileName) = @_;
  print( "Reading features from $fileName...\n");
  # counter statement
  my $countFeatures = 0; my $mRNAcount=0; my $tRNAcount=0; my $exonCount=0; my %CDScount; my %UTRcount; my $nbExonExpanded=0;
  my %UTR3count; my $sizeUTR3=0; my %UTR5count; my $sizeUTR5=0; my %UTRbothSideCount; my %UTR3SideCount; my %UTR5SideCount; my %UTRanyideCount;
  # hash for duplication check statement
  my %geneDupli; my %mRNAdupli; my %tRNAdupli; my %CDSdupli; my %exonDupli; my %UTRdupli;  my %UTR3dupli;  my %UTR5dupli;
  # repeat variables statement
  my $sizeRepeatMasker=0; my %RepeatMaskerDupli; my $sizeRepeatRunner=0; my %RepeatRunnerDupli; my %RepeatMatchPartDupli;
  # hash of data (from feature) statement
  my %genes; my %mRNAHash; my %exonHash; my %cdsHash; my %UTR5Hash; my %UTR3Hash; my %UTRHash; my %tRNAHash; my %repeatHash; my %repeatMatch_part; my %pieceStudied;
  # various variable statement
  my $duplicate="no"; my $exonExpandable=0; my $exonBelongsToMultipleParent=0;

  # read file and decompose it
  while (my $feature = $file_in->read_feature() ) {
      $countFeatures++;
      my $type = $feature->feature_type();
      my $ID = $feature->get_attribute('ID');
      my $seqname=$feature->seqname();
      my $start=$feature->start();
      my $end=$feature->end();
      my $strand=$feature->strand();      

      ####################################################
      ########## Manage feature WITHOUT parent ###########
      ####################################################
      if(! ($feature->has_parent())){
        print "$type $feature\n";
        if ( lc($type) eq 'gene' ) {
          $geneDupli{$ID}++;
          if($geneDupli{$ID} == 1){
            $genes{$ID} = $feature;
            my $pieceStudiedName="$seqname.$strand";
            push( @{ $pieceStudied{$pieceStudiedName}}, $ID );
          }
        next();
        }

        if ( (lc ($feature->source()) =~ m/repeat/) and (lc ($type) eq "match")){
          $RepeatMaskerDupli{$ID}++;
          if ( $RepeatMaskerDupli{$ID} == 1) {
              $repeatHash{$ID} = $feature;
              $sizeRepeatMasker=$sizeRepeatMasker+($feature->end() - $feature->start());
          }
        next();
        }
        if ( (lc ($feature->source()) =~ m/repeatrunner/) and (lc ($type) eq "protein_match")){
          $RepeatRunnerDupli{$ID}++;
          if($RepeatRunnerDupli{$ID} == 1){
              $repeatHash{$ID} = $feature;
              $sizeRepeatRunner=$sizeRepeatRunner+($feature->end() - $feature->start());
          }
        next();
        }
        printf( STDERR "Feature $type not yet taken in account...\n");
      }


      ################################################
      ########## Manage feature WITH parent ##########
      ################################################

      else { # IF NOT a GENE FEATURE
        my $parentID; my @Parent;
        # Manage Attributes
        my $uniqParent="yes";
        my $attributes=$feature->attributes();
        if (ref($feature->get_attribute('Parent')) eq 'ARRAY'){
          @Parent = @{$feature->get_attribute('Parent')};
          $parentID=$Parent[0]; #If several parent ID We take only the first one !
          $uniqParent="no";
          if(lc($type) ne "exon") { # check if multiple feature other than exon. It is not till implemented 
            print "STOP - Your file contains feature $type with multiple parents. Sorry but currently the script manage multiple parent only for exon !\n";
            print "This Warning means you cannot use the \"expand\" option. Consequently you cannot use the \"id\" option that use also the expand option.\n"; exit; }
        }
        else {$parentID= $feature->get_attribute('Parent');}

        # all the following attributes must have parentID
        if ($parentID eq ""){
          print "No Parent attributes found for $ID ! It is mandatory...";exit;
        }

        if ((lc($type) eq 'mrna') || (lc($type) eq 'transcript')){
          if ( !defined( $genes{$parentID}) ) {
            printf( STDERR "ID ".$parentID." dont exists. Gene should be read before ...\n"); exit();
          }
          $mRNAcount++;$mRNAdupli{$ID}++;
          if($mRNAdupli{$ID} == 1){
            addDataToHashOfHash(\%mRNAHash,$parentID,$feature);
          }
          next();
        }
        elsif(lc ($type) eq "trna"){
          $tRNAcount++;$tRNAdupli{$ID}++;
          if($tRNAdupli{$ID} == 1){
            addDataToHashOfHash(\%tRNAHash,$parentID,$feature);
          }
          next();
        }
        elsif (lc($type) eq "exon"){
          if ($uniqParent eq "no"){ #<<<<<<<<<<<<<<<< Expand exon >>>>>>>>>>>>>>>     
            $exonCount++;$exonDupli{$ID}++;
            if($exonDupli{$ID} == 1){
              # for each parent of the exon, create an exon (even if only one...)
              foreach my $mRNA_ID (@Parent){
                $nbExonExpanded++;
                my $featureSaved=$feature->copy();
                
                # Remove all parents with other IDs than $mRNA_ID from $feature.
                foreach my $parent_to_remove ( @{ $featureSaved->parents() } ) {
                  if ( $parent_to_remove->get_attribute('ID') ne $mRNA_ID ) {
                      $featureSaved->remove_parent($parent_to_remove);
                   }
                }

                #change ID to be Uniq
                $featureSaved->set_attribute('ID',"$ID-$mRNA_ID");
                addDataToHashOfHash(\%exonHash,$mRNA_ID,$featureSaved);                       
              }
              $nbExonExpanded--;
            }
          } # end option expand
          else{ # case where no expand option
            if($uniqParent eq "no"){ $exonExpandable=$exonExpandable+$#Parent; $exonBelongsToMultipleParent++; }
            $exonCount++;$exonDupli{$ID}++;
            if($exonDupli{$ID} == 1){
              addDataToHashOfHash(\%exonHash,$parentID,$feature);
            }
          }
        next();
        }
        elsif (lc($type) eq "cds"){ # /!\ CDS feature is described by several features that have the same ID
          my $createdID="$seqname.$start.$end.$ID";
          $CDScount{$ID}++;$CDSdupli{$createdID}++;
          if($CDSdupli{$createdID} == 1){
            addDataToHashOfHash(\%cdsHash,$parentID,$feature);
          }
          next();
        }
        elsif (lc($type) eq "five_prime_utr"){ # /!\ CDS feature is described by several features that have the same ID
          my $createdID="$seqname.$start.$end.$ID";
          $UTR5count{$parentID}++;$UTR5dupli{$createdID}++;
          if($UTR5dupli{$createdID} == 1){
            addDataToHashOfHash(\%UTR5Hash,$parentID,$feature);
            $sizeUTR5=$sizeUTR5+($feature->end() - $feature->start());  
            my $mRNAfeature=$feature->parent();
            my $geneID=$mRNAfeature->get_attribute('Parent');
            $UTR5SideCount{$geneID}++;$UTRanyideCount{$geneID}++;
            if (exists ($UTR3count{$parentID})){
              $UTRbothSideCount{$geneID}++;
            }
          } 
          next();
        }
        elsif (lc($type) eq "three_prime_utr"){ # /!\ CDS feature is described by several features that have the same ID
          my $createdID="$seqname.$start.$end.$ID";
          $UTR3count{$parentID}++;$UTR3dupli{$createdID}++;
          if($UTR3dupli{$createdID} == 1){
            addDataToHashOfHash(\%UTR3Hash,$parentID,$feature); 
            $sizeUTR3=$sizeUTR3+($feature->end() - $feature->start());
            my $mRNAfeature=$feature->parent();
            my $geneID=$mRNAfeature->get_attribute('Parent');
            $UTR3SideCount{$geneID}++;$UTRanyideCount{$geneID}++;
            if (exists ($UTR5count{$parentID})){
              $UTRbothSideCount{$geneID}++;
            }
          } 
          next();
        }
        elsif (lc($type) eq "utr"){ # /!\ CDS feature is described by several features that have the same ID
          my $createdID="$seqname.$start.$end.$ID";
          $UTRcount{$parentID}++;$UTRdupli{$createdID}++;
          if($UTRdupli{$createdID} == 1){
            addDataToHashOfHash(\%UTRHash,$parentID,$feature);
          }      
          next();
        }   
        elsif( (lc ($feature->source()) =~ m/repeat/) and (lc ($type) eq "match_part")){
          $RepeatMatchPartDupli{$ID}++;
          if($RepeatMatchPartDupli{$ID} == 1){
            addDataToHashOfHash(\%repeatMatch_part,$parentID,$feature);
          } 
        next();
        }
        else{
            printf( STDERR "Feature $type not yet taken in account...\n");
#            exit();
        }
      }
  }
  
  $ref_istream->close();

  # Display information for the user:
  my $stringPrint;
  $stringPrint  = "Read $countFeatures features.\n\n" ;
  $stringPrint .= "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n";
  $stringPrint .= "vvvvvvvv Checking of duplicated features vvvvvvvv\n\n";
# About duplicate #
  my $nbDupli=0;
  foreach my $id ( keys %geneDupli ){    
    if ( $geneDupli{$id} != "1"){
      $nbDupli=$nbDupli+($geneDupli{$id}-1);$duplicate="yes";
    }
  }
  $stringPrint .= "Gene => $nbDupli duplicated (ID analyzed)\n";
  $nbDupli=0;
  foreach my $id ( keys %mRNAdupli ){    
    if ( $mRNAdupli{$id} != "1"){
      $nbDupli=$nbDupli+($mRNAdupli{$id}-1);$duplicate="yes";
    }
  }
  $stringPrint .= "mRNA => $nbDupli duplicated (ID analyzed)\n";
  $nbDupli=0;
  foreach my $id ( keys %tRNAdupli ){    
    if ( $tRNAdupli{$id} != "1"){
      $nbDupli=$nbDupli+($tRNAdupli{$id}-1);$duplicate="yes";
    }
  }
  $stringPrint .= "tRNA => $nbDupli duplicated (ID analyzed)\n";
  $nbDupli=0;
  foreach my $id ( keys %exonDupli ){    
    if ( $exonDupli{$id} != "1"){
      $nbDupli=$nbDupli+($exonDupli{$id}-1);$duplicate="yes";
    }
  }
  $stringPrint .= "Exon => $nbDupli duplicated (ID analyzed)\n";
  $nbDupli=0;
  foreach my $id ( keys %CDSdupli ){    
    if ( $CDSdupli{$id} != "1"){
      $nbDupli=$nbDupli+($CDSdupli{$id}-1);$duplicate="yes";
    }
  }
  $stringPrint .= "CDS  => $nbDupli duplicated (seqname().start().end().ID analyzed)\n";
  $nbDupli=0;
  foreach my $id ( keys %UTRdupli ){    
    if ( $UTRdupli{$id} != "1"){
      $nbDupli=$nbDupli+($UTRdupli{$id}-1);$duplicate="yes";
    }
  }
  $stringPrint .= "UTR  => $nbDupli duplicated (seqname().start().end().ID analyzed)\n";
  $nbDupli=0;
  foreach my $id ( keys %UTR3dupli ){    
    if ( $UTR3dupli{$id} != "1"){
      $nbDupli=$nbDupli+($UTR3dupli{$id}-1);$duplicate="yes";
    }
  }
  $stringPrint .= "UTR3 => $nbDupli duplicated (seqname().start().end().ID analyzed)\n";
  $nbDupli=0;
  foreach my $id ( keys %UTR5dupli ){    
    if ( $UTR5dupli{$id} != "1"){
      $nbDupli=$nbDupli+($UTR5dupli{$id}-1);$duplicate="yes";
    }
  }
  $stringPrint .= "UTR5 => $nbDupli duplicated (seqname().start().end().ID analyzed)\n";
  $nbDupli=0;
  foreach my $id ( keys %RepeatMaskerDupli ){    
    if ( $RepeatMaskerDupli{$id} != "1"){
      $nbDupli=$nbDupli+($RepeatMaskerDupli{$id}-1);$duplicate="yes";
    }
  }
  $stringPrint .= "RepeatMasker => $nbDupli duplicated (ID analyzed)\n";
  $nbDupli=0;
  foreach my $id ( keys %RepeatRunnerDupli ){    
    if ( $RepeatRunnerDupli{$id} != "1"){
      $nbDupli=$nbDupli+($RepeatRunnerDupli{$id}-1);$duplicate="yes";
    }
  }
  $stringPrint .= "RepeatRunner => $nbDupli duplicated (ID analyzed)\n";
  $nbDupli=0;
  foreach my $id ( keys %RepeatMatchPartDupli ){    
    if ( $RepeatMatchPartDupli{$id} != "1"){
      $nbDupli=$nbDupli+($RepeatMatchPartDupli{$id}-1);$duplicate="yes";
    }
  }
  $stringPrint .= "RepeatMatchpPart => $nbDupli duplicated (ID analyzed)\n";
  if ($duplicate eq "yes"){ 
    $stringPrint .= "##################################\n# Achthung /\\ Attention /\\ Be carefull => ID duplicate found ! #\n".
      "# Duplicated features have been removed (Keep only one per ID)\n##################################\n\n"; 
  }
  else {$stringPrint .= "##################################\n# Congratulation no duplicated ID #\n##################################\n";}
  $stringPrint .= "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n\n";

# STATISTICS
$stringPrint .=  "vvvvvvvvvvvvvvvvvvvvvvvvvvvv\n";
$stringPrint .=  "vvvvvvvv STATISTICS vvvvvvvv\n\n";
# GENE #
  $stringPrint .= "Gene information:\n";
  my $countGeneUniq = keys (%geneDupli);
  $stringPrint .=  "Read $countGeneUniq total uniq gene\n" ;
  my $nbrmRNAkey = keys (%mRNAHash);
  $stringPrint .= "read gene from mRNA: $nbrmRNAkey\n" ;
  my $nbtRNAkey = keys (%tRNAHash);
  $stringPrint .= "read gene from tRNA: $nbtRNAkey\n\n" ;
# mRNA #
  $stringPrint .= "mRNA information:\n";
  my $nbtotalRNA= ($tRNAcount+$mRNAcount);
  my $nbmRNAisoform=($mRNAcount+$tRNAcount)-$countGeneUniq;
  $stringPrint .= "Read $mRNAcount mRNA. \n";
  $stringPrint .= "Read $tRNAcount tRNA.\n";
  $stringPrint .= "Total RNA= $nbtotalRNA rna\n$nbmRNAisoform mRNA isoforms.\n";
  my $nbrExonKey = keys (%exonHash);
  my $nbRNAwithoutExonDueToMutipleParentalExon=$nbtotalRNA-$nbrExonKey;
  my $verifSize=$nbRNAwithoutExonDueToMutipleParentalExon+$nbrExonKey;
  $stringPrint .= "Nb mRNA key from exons: $nbrExonKey. $nbRNAwithoutExonDueToMutipleParentalExon mRNA without exon du to multiple parents. Total is $verifSize. (Must be the same value as Total RNA)\n\n";  
# EXON #
  $stringPrint .= "exon information:\nRead $exonCount exon \n";
  $stringPrint .= "\n";
# CDS #
  my $nbCDS = keys (%CDScount); 
  $stringPrint .= "CDS information:\nRead $nbCDS CDS\n" ;
  my $nbCDSkey = keys (%CDSdupli); 
  $stringPrint .= "Read $nbCDSkey CDS exon\n\n" ;
# UTR #
  my $nbUTR5exon = keys (%UTR5dupli); my $nbUTR5 = keys (%UTR5count);
  my $nbUTR3exon = keys (%UTR3dupli); my $nbUTR3 = keys (%UTR3count);
  my $nbUTRexon = keys (%UTRdupli);   my $nbUTR = keys (%UTRcount);
  my $totalUTR=$nbUTR5+$nbUTR3+$nbUTR; 
  my $UTRbothSideNb= keys (%UTRbothSideCount);
  my $UTR3SideNB= keys %UTR3SideCount;
  my $UTR5SideNB= keys %UTR5SideCount;
  my $UTRanySideNB= keys %UTRanyideCount;
  $stringPrint .= "UTR information:\nRead $totalUTR UTR exon\n...Read $nbUTR3 UTR3 <=> $nbUTR3exon exon <=> $sizeUTR3 bp length\n...Read $nbUTR5 UTR5 <=> $nbUTR5exon exon <=> $sizeUTR5 bp length\n".
  "...Read $nbUTR UTR <=> $nbUTRexon exon  (whithout more details if come from 3 or 5 prime)\n".
  "...Nb gene that have both UTR: $UTRbothSideNb\n...Nb gene that have 3' UTR: $UTR3SideNB\n...Nb gene that have 5' UTR: $UTR5SideNB\n...Nb gene that have at least one UTR: $UTRanySideNB\n\n" ;
# Repeat #
  $stringPrint .= "Repeat information:\n";
  my $nbRepeatMaskerFeatureUniq = keys (%RepeatMaskerDupli);
  my $nbRepeatRunnerFeatureUniq = keys (%RepeatRunnerDupli);
  my $nbRepeatFeat = $nbRepeatRunnerFeatureUniq+$nbRepeatMaskerFeatureUniq;
  my $nbRepeatMatch = keys (%RepeatMatchPartDupli);
  $stringPrint .= "Read $nbRepeatFeat repeat and $nbRepeatMatch match_part\nWe have $nbRepeatMaskerFeatureUniq repeat features from repeatmasker corresponding to $sizeRepeatMasker bases\n";
  $stringPrint .= "We have $nbRepeatRunnerFeatureUniq repeat features from repeatrunner corresponding to $sizeRepeatRunner bases\n\n";
  $stringPrint .= "Exons expansion added $nbExonExpanded new exons\n\n";
  
  # display
  print "$stringPrint";

  return \%genes,\%mRNAHash,\%exonHash,\%cdsHash,\%UTR5Hash,\%UTR3Hash, \%UTRHash, \%tRNAHash, \%repeatHash, \%repeatMatch_part, \%pieceStudied;
}



__END__

