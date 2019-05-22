#!/usr/bin/env perl

################################################
# maker_checkFusionSplitBetweenTwoBuilds.pl v1 #
# Jacques Dainat 10/2014                       #
# Jacques.dainat@bils.se                       #
################################################

use strict;
use warnings;
use Data::Dumper;
use Carp;
use Getopt::Long;
use IO::File;
use Pod::Usage;

use FindBin qw( $Bin );
use lib "$Bin/../../lib/perl";

use Bio::Tools::GFF;

my $opt_expand="no";my $nbexpand;
my $opt_reffile;
my $opt_tarfile;
my $opt_dirRes;
my $opt_help = 0;

if ( !GetOptions( 'f|ref|reffile=s' => \$opt_reffile,
                  't|tar|tarfile=s' => \$opt_tarfile,
                  'o|out|output=s' => \$opt_dirRes,
                  'h|help!'         => \$opt_help ) )
{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

if ($opt_help) {
    pod2usage( { -verbose => 2,
                 -exitval => 0 } );
}

if ( !( defined($opt_reffile) && defined($opt_tarfile)  && defined($opt_dirRes) ) ) {
    pod2usage( {
           -message => "Must specify 3 parameters:\nReference data gff3 file (--ref) " .
             "\nTargeted gff3 file (--tar)\nOuput directory (--out)",
           -verbose => 0,
           -exitval => 2 } );
}

#####################
# Manage Input File #
#####################
my $ref_istream = IO::File->new();
my $add_istream = IO::File->new();

$ref_istream->open( $opt_reffile, 'r' ) or
  croak(
     sprintf( "Can not open '%s' for reading: %s", $opt_reffile, $! ) );
$add_istream->open( $opt_tarfile, 'r' ) or
  croak(
     sprintf( "Can not open '%s' for reading: %s", $opt_tarfile, $! ) );

my $ref_in = Bio::Tools::GFF->new( -fh => $ref_istream , -gff_version => 3);
my $add_in = Bio::Tools::GFF->new( -fh => $add_istream , -gff_version => 3);
#########################
# END Manage Input File #
#########################
#################################
# Manage Ouput Directory / File #
#################################
if (-d $opt_dirRes){
  print "The output directory choosen already exists. Please give me another Name.\n";exit();
}
my $outDir="";
if ($opt_dirRes =~ /$\//){
  $outDir=$opt_dirRes;
}else{$outDir="$opt_dirRes/";}
unless(mkdir $outDir) {
  die "Unable to create $outDir";
}
my $ostreamResume     = IO::File->new();
my $ostream     = IO::File->new();
my $ostream2    = IO::File->new();
my $ostream3    = IO::File->new();
my $ostream4    = IO::File->new();
my $ostream5    = IO::File->new();
my $ostream6    = IO::File->new();
my $ostream7    = IO::File->new();
my $ostream8    = IO::File->new();

my $opt_outputCluster=$outDir."ouputClusterRef.gff";
my $opt_outputMerge=$outDir."ouputSplitMergeRef.gff";
my $opt_outputCluster2=$outDir."ouputClusterTar.gff";
my $opt_outputMerge2=$outDir."ouputSplitMergeTar.gff";
my $opt_ouputGeneTarMergedInRef=$outDir."ouputGeneTarMergedInRef.gff";
my $opt_ouputGeneRefMergedInTar=$outDir."ouputGeneRefMergedInTar.gff";
my $opt_ouputGeneTarSplitInRef=$outDir."ouputGeneTarSplitInRef.gff";
my $opt_ouputGeneRefSplitInTar=$outDir."ouputGeneRefSplitInTar.gff";
my $resumeFile=$outDir."resume.txt";
$ostreamResume->open( $resumeFile, 'w' ) or
  croak(
      sprintf( "Can not open '%s' for writing %s", $resumeFile, $! )
  );
$ostream->open( $opt_outputCluster, 'w' ) or
  croak(
      sprintf( "Can not open '%s' for writing %s", $opt_outputCluster, $! )
  );
$ostream2->open( $opt_outputMerge, 'w' ) or
  croak(
      sprintf( "Can not open '%s' for writing %s", $opt_outputMerge, $! )
  );
$ostream3->open( $opt_outputCluster2, 'w' ) or
  croak(
      sprintf( "Can not open '%s' for writing %s", $opt_outputCluster2, $! )
  );
$ostream4->open( $opt_outputMerge2, 'w' ) or
  croak(
      sprintf( "Can not open '%s' for writing %s", $opt_outputMerge2, $! )
  );
$ostream5->open( $opt_ouputGeneTarMergedInRef, 'w' ) or
  croak(
      sprintf( "Can not open '%s' for writing %s", $opt_ouputGeneTarMergedInRef, $! )
  );
$ostream6->open( $opt_ouputGeneRefMergedInTar, 'w' ) or
  croak(
      sprintf( "Can not open '%s' for writing %s", $opt_ouputGeneRefMergedInTar, $! )
  );
$ostream7->open( $opt_ouputGeneTarSplitInRef, 'w' ) or
  croak(
      sprintf( "Can not open '%s' for writing %s", $opt_ouputGeneTarSplitInRef, $! )
  );
$ostream8->open( $opt_ouputGeneRefSplitInTar, 'w' ) or
  croak(
      sprintf( "Can not open '%s' for writing %s", $opt_ouputGeneRefSplitInTar, $! )
  );

my $outputCluster = $ostream ;
my $outputMerge = $ostream2 ;
my $outputCluster2 = $ostream3 ;
my $outputMerge2 = $ostream4 ;
my $ouputGeneTarMergedInRef = $ostream5 ;
my $ouputGeneRefMergedInTar = $ostream6 ;
my $ouputGeneTarSplitInRef = $ostream7 ;
my $ouputGeneRefSplitInTar = $ostream8 ;

print $outputCluster "##gff-version 3\n";
print $outputMerge "##gff-version 3\n";
print $outputCluster2 "##gff-version 3\n";
print $outputMerge2 "##gff-version 3\n";
print $ouputGeneTarMergedInRef "##gff-version 3\n";
print $ouputGeneRefMergedInTar "##gff-version 3\n";
print $ouputGeneTarSplitInRef "##gff-version 3\n";
print $ouputGeneRefSplitInTar "##gff-version 3\n";
#####################################
# END Manage Ouput Directory / File #
#####################################

# declaration hashes reference
# data from file one
my $ref_genes; my $ref_trnagenes; my $refmRNA;my $refexon; my $refcds; my $refUTR5; my $refUTR3; my $reftRNA; my $refpieceStudied;
# data from file 2
my $tar_genes; my $tar_trnagenes; my $tarmRNA;my $tarexon; my $tarcds; my $tarUTR5; my $tarUTR3; my $tartRNA; my $tarpieceStudied;
# data created by mixing data from file1 and file2
my $mix_genes; my $mix_trnagenes; my $mixmRNA;my $mixexon; my $mixcds; my $mixUTR5; my $mixUTR3; my $mixtRNA; my $mixpieceStudied;

#Parse GFF genes,\%mRNA,\%exon,\%cds,\%UTR5,\%UTR3,\%tRNA;
($ref_genes,$ref_trnagenes,$refmRNA,$refexon, $refcds, $refUTR5, $refUTR3, $reftRNA, $refpieceStudied) = parseGFF($ref_in,$opt_reffile, $opt_expand);
($tar_genes,$tar_trnagenes,$tarmRNA,$tarexon, $tarcds, $tarUTR5, $tarUTR3, $tartRNA, $tarpieceStudied) = parseGFF($add_in,$opt_tarfile, $opt_expand);
print("Parsing Finished\n");

sortByPos($refexon);
sortByPos($tarexon);
sortByPos($refcds);
sortByPos($tarcds);
sortByPos($refUTR5);
sortByPos($tarUTR5);
sortByPos($refUTR3);
sortByPos($tarUTR3);
print("Sort Finished\n");

my $countGeneA=0;
my $countGeneB=0;
my $nbGeneFromContigSpecificRef=0;
my $nbSpecificToTar=0;
my %refpieceStudiedHash=%$refpieceStudied;
my %tarpieceStudiedHash=%$tarpieceStudied;
my @listContigsBoth;my @listContigsA;my @listContigsB;


                                      ##############################################################
                                      # A) Manage scaffold (stranded) annotated only by one method # print genes and remove them from hashes
                                      ##############################################################

#########################################################
# A.1) Manage gene from contig annotated only in file 1 # contig direction took in account
#########################################################
my $method1=0;
foreach my $ContigName (keys %refpieceStudiedHash) {
  if (! exists $tarpieceStudiedHash{$ContigName}){
    $method1++;
    my $gene=$refpieceStudiedHash{$ContigName};
    foreach my $geneName (@{$gene}) {
      $nbGeneFromContigSpecificRef++;
      $countGeneA++;
    }
    my $tmpContigName = $ContigName;
    $tmpContigName =~ s/[+-]//g;
    push(@listContigsA, $tmpContigName);
    delete $refpieceStudiedHash{$ContigName};
  }
}
my @listContigsAUniq = Array_Unique(@listContigsA);
my $nbContigAUniq=$#listContigsAUniq+1;
print("$nbContigAUniq contigs annotated only in the reference build\n");
#^^^^^^^^^^^^^^^^^^^  END ^^^^^^^^^^^^^^^^^^^

#########################################################
# A.2) Manage gene from contig annotated only in file 2 # contig direction took in account
#########################################################
my $method2=0;
foreach my $ContigName (keys %tarpieceStudiedHash) {
  if (! exists $refpieceStudiedHash{$ContigName}){
    $method2++;
    my $gene=$tarpieceStudiedHash{$ContigName};
    foreach my $geneName (@{$gene}) {
      $nbSpecificToTar++;
      $countGeneB++;
    }
    my $tmpContigName = $ContigName;
    $tmpContigName =~ s/[+-]//g;
    push(@listContigsB, $tmpContigName);
    delete $tarpieceStudiedHash{$ContigName};
  }
}
my @listContigsBUniq = Array_Unique(@listContigsB);
my $nbContigBUniq=$#listContigsBUniq+1;
print("$nbContigBUniq contigs annotated only in the target build\n");
#^^^^^^^^^^^^^^^^^^^  END ^^^^^^^^^^^^^^^^^^^

                                      ############################################################################
                                      # B) Manage scaffold (stranded) which contain annotation from both methods # 
                                      ############################################################################
my $countSingleGeneB=0;
my $OverlapingB=0;
my $OverlapingA=0;
my $nbNoOverlapingA=0;
my $nbNoOverlapingB=0;
my %toStretch;
my %clusterCase;
my %fusionORsplit;
my $fusion;
my $split;
print "Now we are analyzing the contigs containing annotations from both builds:\n";
foreach my $ContigName (keys %refpieceStudiedHash) {
#  print("Study of ContigName $ContigName\n");
  my $tmpContigName = $ContigName;
  $tmpContigName =~ s/[+-]//g;
  push(@listContigsBoth, $tmpContigName);
  #for each gene A to B
  my @tempTabName = @{$refpieceStudiedHash{$ContigName}}; # use of temporary variable to be sure to loop over all element.
  foreach my $GeneName (@tempTabName) { 
    $countGeneA++;
    #print("\nStudy overlap of gene $GeneName\n");
   
    ###### Test if gene already studied (started by another gene but due to overlap it is already studied)
    my $geneAlreadyStudied="yes";
    foreach my $gene (@{$refpieceStudiedHash{$ContigName}}){
      if ($gene eq $GeneName){$geneAlreadyStudied="no";last;} 
    }
    if ($geneAlreadyStudied eq "yes"){next;}#    print "$GeneName already analyzed because overlap other chunck: NEXT  \n";
    ##### End test if already studied. If not, we continue
   
    else
    {
      # Declare table which I will work with
      my @ListRefOverlapAtotest;my @ListOverlapAtested;my @ListNoOverlapA; my @ListOverlapAtestnoneed; my @ListPerfectOverlapA;
      my @ListRefOverlapBtotest;my @ListOverlapBtested;my @ListNoOverlapB; my @ListOverlapBtestnoneed; my @ListPerfectOverlapB;

    my @LinkTocurrentGeneFeature=$ref_genes->{$GeneName};  #### >>>>>> BASE OF THE FEATURE TESTED FOR OVERLAP !!! CURRENTLY WE CHECK THE GENE FEATURE !
    my $firstRound=0;

    push (@ListRefOverlapAtotest, @LinkTocurrentGeneFeature);

      while ($#ListRefOverlapAtotest != -1 or $#ListRefOverlapBtotest != -1){

            $firstRound++;
            if ($#ListRefOverlapAtotest != -1){
              # test every A
              my ($ListRefOverlapBtotest, $ListOverlapBtestnoneed, $ListOverlapAtested, $ListNoOverlapA, $ListPerfectOverlapA) = retrieveAllOverlap( $firstRound, $refpieceStudiedHash{$ContigName}, $tarpieceStudiedHash{$ContigName}, $tar_genes, \@ListRefOverlapAtotest, \@ListOverlapAtested, \@ListNoOverlapA, \@ListPerfectOverlapA, \@ListOverlapBtestnoneed,1);
              @ListRefOverlapBtotest = @$ListRefOverlapBtotest;
      #        print "Nb overlap B to test = $#ListRefOverlapBtotest\n";
              @ListOverlapBtestnoneed = @$ListOverlapBtestnoneed;
              @ListRefOverlapAtotest = ();
     #         print Dumper($ListOverlapAtested);
     #         foreach my $i (@{$ListOverlapAtested->[0]}){print "\nPPPPP $i\n"; }
     #         print "\nPPPPP @{$ListOverlapAtested->[0]}\n";
              @ListOverlapAtested = @$ListOverlapAtested;
              @ListPerfectOverlapA = @$ListPerfectOverlapA;
              @ListNoOverlapA = @$ListNoOverlapA;

      #        print "ListOverlapAtested A $#ListRefOverlapBtotest, $#ListOverlapBtestnoneed, $#ListOverlapAtested, $#ListNoOverlapA, $#ListPerfectOverlapA\n";
              next(); #stop here and avoid test B
            }
           # print "\nPPPPP1 $ListOverlapAtested[0]\n";
            if ($#ListRefOverlapBtotest != -1){
              # test every B
              #print "List size overlap B to test = $#ListRefOverlapBtotest $firstRound \n";

              my ($ListRefOverlapAtotest, $ListOverlapAtestnoneed, $ListOverlapBtested, $ListNoOverlapB, $ListPerfectOverlapB) = retrieveAllOverlap( $firstRound, $tarpieceStudiedHash{$ContigName}, $refpieceStudiedHash{$ContigName}, $ref_genes, \@ListRefOverlapBtotest, \@ListOverlapBtested,\@ListNoOverlapB, \@ListPerfectOverlapB, \@ListOverlapAtestnoneed,0);                        
              @ListRefOverlapAtotest = @$ListRefOverlapAtotest;
              @ListOverlapAtestnoneed = @$ListOverlapAtestnoneed;

              @ListRefOverlapBtotest = ();
              @ListOverlapBtested = @$ListOverlapBtested;
              @ListPerfectOverlapB =  @$ListPerfectOverlapB;           
              @ListNoOverlapB = @$ListNoOverlapB;
       #       print "ListOverlapBtested B $#ListOverlapBtested\n";
            #  print "\nPPPPP2 @{$ListOverlapAtested[0]}\n";
              next(); #stop here
            }
          }
      $nbNoOverlapingA=$nbNoOverlapingA+$#ListNoOverlapA+1;
      $OverlapingA=$OverlapingA+$#ListOverlapAtested+$#ListOverlapAtestnoneed+$#ListPerfectOverlapA+3;

      $nbNoOverlapingB=$nbNoOverlapingB+$#ListNoOverlapB+1;
      $OverlapingB=$OverlapingB+$#ListOverlapBtested+$#ListOverlapBtestnoneed+$#ListPerfectOverlapA+3;

   
      my $nbFragment = ($#ListOverlapAtested+$#ListOverlapBtested+2);
      #print ("\nOVERLAP step 1 END\n");
      if ( $nbFragment < 2 ){
  ##################################
  # Manage single gene from file 1 #
  ##################################
         #HEre can be printed => Single gene from file A and PerfectMatch from A or B depending to an iption like my $contigFusionA = "ok";  < /!\ >
  #^^^^^^^^^^^^  END ^^^^^^^^^^^^
      }
      elsif ( $nbFragment == 2){ #case can be stretched  
      #print "ListOverlapAtested @ListOverlapAtested ListOverlapBtested @ListOverlapBtested";
        push ( @{ $toStretch{$GeneName} }, [@ListOverlapAtested],[@ListOverlapBtested]);
      }
      elsif ($nbFragment == 3){
        #print "FUSION (ref->tar)/SPLIT(tar->ref) case\n";
        push ( @{ $fusionORsplit{$GeneName} }, [[@ListOverlapAtested],[@ListOverlapBtested]]);
        # Fusion in target Build
        if ($#ListOverlapAtested > $#ListOverlapBtested){
          $fusion++;
        } else{$split++;}  # Split in the target Build
      } 
      else{ # Cluster case (more than 3 segments)
  #      print "CLUSTER CASE";
        push ( @{ $clusterCase{$GeneName} }, [[@ListOverlapAtested],[@ListOverlapBtested]]); 
      }

    }
  }

##################################
# Manage single gene from file 2 #
##################################
my @ListSingleB=@{$tarpieceStudiedHash{$ContigName}};
$countSingleGeneB= $countSingleGeneB + $#ListSingleB+1;  
## END ##
}


#################
# Display results

my $nbrSplitMergeToManage = keys (%fusionORsplit);
my $nbrClusterToManage = keys (%clusterCase);
my $nbrStretchingCaseToManage = keys (%toStretch);

my $totalSpecificToA=$nbGeneFromContigSpecificRef+$nbNoOverlapingA;
my $totalA=$OverlapingA+$totalSpecificToA;

my $totalSpecificToB=$nbSpecificToTar+$nbNoOverlapingB+$countSingleGeneB;
my $totalB=$OverlapingB+$totalSpecificToB;

my @listContigsUniq = Array_Unique(@listContigsBoth);
my $nbContigUniq=$#listContigsUniq+1;

my $resultToPrint="";
$resultToPrint.= "\n\n######### RESULTS #########:\n\n";
$resultToPrint.=  "File1 ($opt_reffile):\n";
$resultToPrint.=  "We studied ($nbContigUniq) contigs and $countGeneA genes\n";
$resultToPrint.=  "=> $totalSpecificToA genes are nonoverlapping (i.e unique to this gene build => $nbGeneFromContigSpecificRef are on contig(s) only annotated in this gene build and $nbNoOverlapingA are on contigs also annotated by the second gene build.)\n";
$resultToPrint.=  "=> $OverlapingA genes overlap genes from file 2\n";
$resultToPrint.=  "Total gene checked : $totalA \n\n";
$resultToPrint.=  "File2 ($opt_tarfile):\n";
$resultToPrint.=  "=> $totalSpecificToB genes are nonoverlapping (i.e unique to this gene build => $nbSpecificToTar are on contig(s) only annotated in this gene build and $countSingleGeneB are on contigs also annotated by the first gene build.)\n";
$resultToPrint.=  "=> $OverlapingB genes overlap genes from file 1\n";
$resultToPrint.=  "Total gene checked = $totalB\n\n";
$resultToPrint.=  "Results:\n";
$resultToPrint.=  "Number of SPLIT/MERGE case : $nbrSplitMergeToManage\n";
my $nbGeneImplicatedF=$fusion*2;
my $nbGeneImplicatedS=$split*2;
$resultToPrint.=  "More precisely: $fusion cases of fusion in the target build detected. (corresponding to $nbGeneImplicatedF genes from reference build implicated) \n";
$resultToPrint.=  "                $split cases of split in the target build detected. (corresponding to $nbGeneImplicatedS genes from target build implicated) \n";
$resultToPrint.=  "Number of CLUSTER case : $nbrClusterToManage\n";
$resultToPrint.=  "=> Result are written in gff3 format in $outDir directory\n\n";
print $resultToPrint;
print $ostreamResume "$resultToPrint";

################################
# B.1) Manage case to stretch  #
################################

# Nothing to do

##################################
# B.2) Manage split/fusion case  #
##################################
for my $geneKey (keys %fusionORsplit){
  my @A_geneList = @{ @{ $ { $fusionORsplit{$geneKey} } [0] } [0]};
  my @B_geneList = @{ @{ $ { $fusionORsplit{$geneKey} } [0] } [1]};
   # print all gene A
  for my $geneInfo (@A_geneList){
    my $geneName = @{$geneInfo}[0];
    printgenes($outputMerge, $geneName, $ref_genes, $refmRNA, $refexon, $refcds, $refUTR5, $refUTR3, $reftRNA) ;
  }
  # print all gene B
  for my $geneInfo (@B_geneList){
    my $geneName = @{$geneInfo}[0];
    printgenes($outputMerge2, $geneName, $tar_genes, $tarmRNA, $tarexon, $tarcds, $tarUTR5, $tarUTR3, $tartRNA) ;
  }
  #Print Independently Fusion or Split
  if ($#A_geneList > $#B_geneList){
    #Print ref genes merged in target build
    for my $geneInfo (@A_geneList){
      my $geneName = @{$geneInfo}[0];
      printgenes($ouputGeneRefMergedInTar, $geneName, $ref_genes, $refmRNA, $refexon, $refcds, $refUTR5, $refUTR3, $reftRNA) ;
    }
    # print result of gene fusion in target build
    for my $geneInfo (@B_geneList){
      my $geneName = @{$geneInfo}[0];
      printgenes($ouputGeneTarSplitInRef, $geneName, $tar_genes, $tarmRNA, $tarexon, $tarcds, $tarUTR5, $tarUTR3, $tartRNA) ;
    }
  }
  else{
    #Print target genes merged in ref build
    for my $geneInfo (@A_geneList){
      my $geneName = @{$geneInfo}[0];
      printgenes($ouputGeneTarMergedInRef, $geneName, $ref_genes, $refmRNA, $refexon, $refcds, $refUTR5, $refUTR3, $reftRNA) ;
    }
    # print result of gene fusion in ref build 
    for my $geneInfo (@B_geneList){
      my $geneName = @{$geneInfo}[0];
      printgenes($ouputGeneRefSplitInTar, $geneName, $tar_genes, $tarmRNA, $tarexon, $tarcds, $tarUTR5, $tarUTR3, $tartRNA) ;
    }
  }
}

#############################
# B.2) Manage Cluster case  #
#############################
for my $geneKey (keys %clusterCase){
#  print "/!\\WARNING/!\\ Cluster case ! We keep the genes intact but you have to manage manualy this case:\n";
#  print "$geneKey\n";
  my @A_geneList = @{ @{ $ { $clusterCase{$geneKey} } [0] } [0]};
  my @B_geneList = @{ @{ $ { $clusterCase{$geneKey} } [0] } [1]};
  # print all gene A
  for my $geneInfo (@A_geneList){
    my $geneName = @{$geneInfo}[0];
    printgenes($outputCluster, $geneName, $ref_genes, $refmRNA, $refexon, $refcds, $refUTR5, $refUTR3, $reftRNA) ;
  }
  # print all gene B
  for my $geneInfo (@B_geneList){
    my $geneName = @{$geneInfo}[0];
    printgenes($outputCluster2, $geneName, $tar_genes, $tarmRNA, $tarexon, $tarcds, $tarUTR5, $tarUTR3, $tartRNA) ;
  } 
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

sub Array_Unique
{
    my @List = @_;
    my %FutureList;
    foreach(@List)
    {
        $FutureList{$_} = 1; # remove doublons 
    }
    return (keys(%FutureList));
}

sub retrieveAllOverlap {
  my ($firstRound, $pieceStudiedHashA, $pieceStudiedHashB, $genesHashB,  $ListA_Overlap_totest, $ListA_Overlap_tostretch, $ListA_NoOverlap, $ListA_PerfectOverlap, $ListOverlaptestnoneedB,$ok) = @_;
  my @newListrefOverlaptotestB;
#  my $sizeListrefOverlaptotestA = @$ListrefOverlaptotestA; my $sizeListOverlaptestedA = @$ListOverlaptestedA; my $sizeListNoOverlapA = @$ListNoOverlapA; my $sizeListOverlaptestnoneedB = @$ListOverlaptestnoneedB; my $sizeListPerfectOverlapA = @$ListPerfectOverlapA;
#  print "firstRound $firstRound\tListrefOverlaptotestA $sizeListrefOverlaptotestA\tListOverlaptestedA $sizeListOverlaptestedA\tListNoOverlapA $sizeListNoOverlapA\tListOverlaptestnoneedB $sizeListOverlaptestnoneedB\tListPerfectOverlapA $sizeListPerfectOverlapA\n";
#  print "Begin retrieveAllOverlap\n";
  for my $linkToGene (@{$ListA_Overlap_totest}) {
      #print "Link = @{$linkToGene} \n";
      my $startA=@{$linkToGene}[4];
      my $endA=@{$linkToGene}[5];
      my $wecanremove="no";my @listGneToRemove;
      my $Overlaped="no";
      my $PerfectOverlap="no";

      for my $geneOppo (@{$pieceStudiedHashB}){ # allow to work only on gene on the contig
        my $startB=$genesHashB->{$geneOppo}[4];
        my $endB=$genesHashB->{$geneOppo}[5];
        
        my $resuOverlap = testOverlap($startA, $endA, $startB, $endB);  ## <====== Call Overlap method
        if ($resuOverlap eq "perfectOverlap"){
          push (@{$ListA_PerfectOverlap}, [@{$linkToGene}] ); # keep reference
          $wecanremove="yes";
          push (@listGneToRemove, $geneOppo);  ################## <=============== Change to keep the possible fusion ??
          $Overlaped="yes";$PerfectOverlap="yes";
          #print "OVERLAP FOUND => perfect Overlap\n";
          last;
        }
        elsif ($resuOverlap eq "noNeedToverify"){
         # print "geneOppo @{$genesHash->{$geneOppo}}\n";
          push (@{$ListOverlaptestnoneedB}, [ @{$genesHashB->{$geneOppo}} ] ); # save features and remove from HashGeneral
          $wecanremove="yes"; # Need to be out of the loop for remove the gene because affect the loop
          push (@listGneToRemove, $geneOppo);   
          $Overlaped="yes";
          #print "OVERLAP FOUND => No NeedToverify\n";
        }
        elsif ($resuOverlap eq "needToVerify"){
         # print "geneOppo $geneOppo\n";
          push (@newListrefOverlaptotestB, $genesHashB->{$geneOppo} );
          $Overlaped="yes";
          #print "OVERLAP FOUND => needToVerify $genesHashB->{$geneOppo}\n";
        }
      }

     if($Overlaped eq "no"){ 
        if ($firstRound == "1"){
          #print "NO Overlap\n";
          push (@{$ListA_NoOverlap}, [@{$linkToGene}]);
        }
        else{ 
          #print "No More Overlap\n";
          push (@{$ListA_Overlap_tostretch}, [@{$linkToGene}]);
        }
      } 
     elsif ($PerfectOverlap ne "yes"){
       #print "test overlap stretched @{$linkToGene}";
        push (@{$ListA_Overlap_tostretch}, [@{$linkToGene}]);
      }

      if($wecanremove eq "yes") {  
        for my $gene (@listGneToRemove){
        removeElementInList($pieceStudiedHashB, $gene);
       }
      }
      # REMOVE
      # Need to be deleted to not re-use it for retrieve overlap => Because we will found again the same // The removed one will be display through @ListA_Overlap_tested
      removeElementInList($pieceStudiedHashA, $linkToGene->[0]);
  }
  return (\@newListrefOverlaptotestB, $ListOverlaptestnoneedB, $ListA_Overlap_tostretch, $ListA_NoOverlap, $ListA_PerfectOverlap);
}

# This method allows to shift on left all number key of a hash
sub deleteFirstElementAndReorganise {
  my ($refHash) = @_;

  my %tmpHash=%$refHash;
  my $lastPos= keys %tmpHash;

  delete $refHash->{$lastPos};
  foreach my $key (keys %{$refHash}){
    $refHash->{$key} = $tmpHash{$key+1};
  }
}

sub firstExonRemovedAffectedCDS{
  my ($exon,$cds)=@_;
  my $exonStart=$exon->[4];
  my $exonEnd=$exon->[5];
  my $cdsStart=$cds->[4];
  my $cdsEnd=$cds->[5];
  my $affected="";
  if (($cdsStart >= $exonStart) and ($cdsEnd >= $exonEnd)) {
    $affected="no";
  }
  else{$affected="yes";}
return $affected;
}

sub lastExonRemovedAffectedCDS{
  my ($exon,$cds)=@_;
  my $exonStart=$exon->[4];
  my $exonEnd=$exon->[5];
  my $cdsStart=$cds->[4];
  my $cdsEnd=$cds->[5];
  my $affected="";
  if (($cdsStart <= $exonStart) and ($cdsEnd <= $exonEnd)) {
    $affected="no";
  }
  else{$affected="yes";}
return $affected;
}

sub exonFinishByCDS{
  my ($exon,$cds)=@_;
  my   $exonEnd=$exon->[5];
  my   $cdsEnd=$cds->[5];
  my $finishByCDS="";
  if ( $exonEnd == $cdsEnd ) {
    $finishByCDS="yes";
  }
  else{$finishByCDS="no";}
return $finishByCDS;
}

sub exonStartByCDS{
  my ($exon,$cds, $strand)=@_;
  my   $exonStart;
  my   $cdsStart;

  if ($strand eq "+"){
     $exonStart=$exon->[4];
     $cdsStart=$cds->[4]; 
  }
  else{
     $exonStart=$exon->[5];
     $cdsStart=$cds->[5];
  }
  my $startByCDS="";
  if ( $exonStart == $cdsStart ) {
    $startByCDS="yes";
  }
  else{$startByCDS="no";}
return $startByCDS;
}

sub featureOnExon{ # start or end should be common !.
  my ($feature1,$feature2)=@_;
  my $feature1Start=$feature1->[4];
  my $feature1End=$feature1->[5];
  my $feature2Start=$feature2->[4];
  my $feature2End=$feature2->[5];
  my $onExon="";
  if (($feature1Start == $feature2Start) or ($feature1End == $feature2End)) {
    $onExon="yes";
  }
  else{$onExon="no";}
return $onExon;
}

sub chooseLongerCDSfrommRNAandPrint{
  my ($output, $A_mRNAname, $geneA, $B_mRNAname, $geneB)= @_;
  my $CDSsizeA=CDSsize($A_mRNAname,"ref");
  my $CDSsizeB=CDSsize($B_mRNAname,"tar");
  if ($CDSsizeA >= $CDSsizeB){
      printgenes($output, $geneA, $ref_genes, $refmRNA, $refexon, $refcds, $refUTR5, $refUTR3, $reftRNA); #<========== # HEre can be printed <== 
  }
  else{
      printgenes($output, $geneB, $tar_genes, $tarmRNA, $tarexon, $tarcds, $tarUTR5, $tarUTR3, $reftRNA)  #<========== # HEre can be printed <== 
  }
}

sub keepmRNAofTheLongestCDS{
  my ($mRNAHash,$type) = @_;
  my $CDSsize=0;
  my $mRNA="";
   foreach my $mRNAlistInfo (keys %{$mRNAHash} ){  # For each mRNA 
    #get mRNA name
    my $mRNAname= $mRNAHash->{$mRNAlistInfo}[0];
    my $mRNAnameSize=CDSsize($mRNAname,$type);
    
    if ($mRNAnameSize > $CDSsize){
      $mRNA=$mRNAname;
    }
  }
  return $mRNA;
}

sub CDSsize {
  my ($mRNAname,$WhichCDS) = @_;

  my $CDS="";my $CDSsize=0;
  if ($WhichCDS eq "ref"){$CDS=$refcds;}else{$CDS=$tarcds;}
  my %hashCDs=@{${$CDS}{$mRNAname}};
  
  foreach my $key (keys %hashCDs){
    my $chunckSize=$hashCDs{$key}[5]-$hashCDs{$key}[4];
    $CDSsize=$CDSsize+$chunckSize;
  }
  return $CDSsize;
}

sub removeElementInList {
  my ($hash1, $element_omitted) = @_;

  @{$hash1}= grep { $_ ne $element_omitted } @{$hash1}; #remove element of the list
#  print "What I delete ?  $hash2->{$element_omitted}";
#  delete $hash2->{$element_omitted};   # remove tuple in hash
#  print "Now I deleted verification:";
#  if (! exists $hash2->{$element_omitted}){
#     print "DELETION OK \n";
#   }
#   else {print "DELETION ERROR \n";}
}

sub testOverlap {
  my ($startA, $endA, $startB, $endB) = @_;
  if($startA == $startB and $endA == $endB){ #overlap perfect     ----
    return "perfectOverlap";                 #                    ----
  }
  elsif($startA <= $startB and $endA >= $endB){# No need to verify  --------
      return "noNeedToverify";                 #                       --
    }
  elsif (($startA >= $startB and $startA <= $endB) or ($endA >= $startB and $endA <= $endB)) {  #   ---      ---   ---
    return "needToVerify";                                                                      # -------  ----     ------
  }
  else{
    return "noOverlap";
  }
}

sub printNewGene{
  my($output, $A_gene, $A_Longest_mRNA, $A_exon, $A_CDS, $A_UTR, $side)= @_;
  print "\nSTART print BRICK\n";
  my $geneName = $A_gene->[0];
  my $geneID = "new_$A_gene->[0]";
  $A_Longest_mRNA->[0] =~ /.*(-mRNA-.*)/ ;
  my $mRNAName = "$geneID$1";

  brickToprintTabNewName($output, $A_gene, $geneID, $geneID);
  brickToprintTabNewName($output, $A_Longest_mRNA, $mRNAName, $geneID);
  brickToprintHashNewName($output, $A_exon, $mRNAName);
  brickToprintHashNewName($output, $A_CDS, $mRNAName);

  if($side eq "right"){ # print UTR5 fron B and UTR3 from A
    if (exists $refUTR5->{$geneName}){ #refUTR is getting as variable reaching fron everywhere
      my %hashUTR5=@{$tarUTR5->{$geneName}};
      brickToprintHashNewName($output, \%hashUTR5, $mRNAName);
    }
    brickToprintHashNewName($output, $A_UTR, $mRNAName);
  }
  elsif($side eq "left"){
    brickToprintHash($output, $A_UTR);
    if (exists $refUTR3->{$geneName}){
      my %hashUTR3=@{$refUTR3->{$geneName}};
      brickToprintHashNewName($output, \%hashUTR3, $mRNAName);
    }
  }
  else{
    #UTR5 and UTR3 data are in the same hash
    brickToprintHashNewName($output, $A_UTR, $mRNAName);
    }
  print "END print BRICK\n\n";
}

sub brickToprintTab{
  my ($output, $refTab)=@_;

  my $cpt=0;
  foreach my $Element (@{$refTab}) {
    if (!($cpt==0 or $cpt==9)){print $output "$Element\t";}
    if ($cpt == 9){ 
      foreach my $element (@{$Element}){
        print $output "$element->[0]=$element->[1];";
      }
      print $output "\n"
    }
    $cpt++;
  }
}

sub brickToprintTabNewName{
  my ($output, $refTab, $ID, $Parent)=@_;

  my $cpt=0;
  foreach my $Element (@{$refTab}) {
    if (!($cpt==0 or $cpt==9)){print $output "$Element\t";}
    if ($cpt == 9){ 
      foreach my $element (@{$Element}){
        if($element->[0] eq "ID"){
          print $output "ID=$ID;";
        }
        elsif($element->[0] eq "Parent"){
          print $output "Parent=$Parent;";
        }
        elsif($element->[0] eq "Name"){
          print $output "Name=$ID;";
        }
        else{
          print $output "$element->[0]=$element->[1];";
        }
      }
      print $output "\n"
    }
    $cpt++;
  }
}

sub brickToprintHashNewName{
  my ($output, $refEle, $prefix)=@_;

  my %refEle = %$refEle;
  my $nbEle = keys %refEle;
  for (my $cptEl = 1 ; $cptEl <= $nbEle ; $cptEl++) { 
    my @tabEle=@{$refEle{$cptEl}};
    my $cpt=0;
  #   # print exon info
    foreach my $Element (@tabEle) {
      if (!($cpt==0 or $cpt==9)){print $output "$Element\t";}
        if ($cpt == 9){
          foreach my $element (@{$Element}){
            if($element->[0] eq "ID"){
            $element->[1] =~ /.+-mRNA-[0-9]*(:.*)/ ;
            print $output "ID=$prefix$1;old_name=$element->[1];";
            }
            elsif($element->[0] eq "Parent"){
              print $output "Parent=$prefix;old_parent=$element->[1];";
            }
            else{
              print $output "$element->[0]=$element->[1];";
            }
          }
      print $output "\n"
      }
    $cpt++;
    }
  }
}

# print in sorted
sub brickToprintHash{
  my ($output, $refEle)=@_;

  my %refEle = %$refEle;
  my $nbEle = keys %refEle;
  for (my $cptEl = 1 ; $cptEl <= $nbEle ; $cptEl++) { 
    my @tabEle=@{$refEle{$cptEl}};
    my $cpt=0;
  #   # print exon info
    foreach my $Element (@tabEle) {
      if (!($cpt==0 or $cpt==9)){print $output "$Element\t";}
      if ($cpt == 9){ 
         foreach my $element (@{$Element}){
            print $output "$element->[0]=$element->[1];";
          }
      print $output "\n"
      }
    $cpt++;
    }
  }
}

#This function use printgene
sub printgenes {
  my ($output, $gene, $OriginalGenesHash,  $mRNA, $exon, $cds, $UTR5, $UTR3, $tRNA) = @_;
  my %OriginalGenesHashOk = %$OriginalGenesHash;
  
  #if $gene is list/hash of gene => print all of them 
  if (ref $gene eq 'ARRAY') {
    foreach my $geneName (@{$gene}) {
      my $current_geneTab = $OriginalGenesHashOk{$geneName};
      printgene($output, $current_geneTab, $mRNA, $exon, $cds, $UTR5, $UTR3, $tRNA);
    }    
  }
  #Print one gene if is just one gene name
  else{
    my $current_geneTab = $OriginalGenesHashOk{$gene};
    printgene($output, $current_geneTab, $mRNA, $exon, $cds, $UTR5, $UTR3, $tRNA);
  }
}

sub printgene {
  my ($output, $genetab, $mRNA, $exon, $cds, $UTR5, $UTR3, $tRNA) = @_;
  my $geneName=$genetab->[0];
  #print gene tab
  brickToprintTab($output, $genetab);
  # Get hash of hash of mRNA 
  if(exists $mRNA->{$geneName}){
    my %mRNAHash=@{$mRNA->{$geneName}};
    # For each mRNA known
    foreach my $mRNAnum (keys %mRNAHash) {
      #Get mRNAhash of the current gene studied
      my $current_mRNAHash=$mRNAHash{$mRNAnum};
      brickToprintTab($output, $current_mRNAHash);
      my $mRNAname=@{$current_mRNAHash}[0];
      #Get exonhash of the current mRNA studied
      if(exists $exon->{$mRNAname}){      # If exon are compacted some mRNA could not have exons ...
        my %exonHash=@{$exon->{$mRNAname}};   
        brickToprintHash($output, \%exonHash);
      }
      #Get cdshash of the current mRNA studied
      my %cdsHash=@{$cds->{$mRNAname}};
      brickToprintHash($output, \%cdsHash);
      #Get UTR5hash of the current mRNA studied if exist
      if (exists $UTR5->{$mRNAname}){
        my %UTR5Hash=@{$UTR5->{$mRNAname}};
        brickToprintHash($output, \%UTR5Hash);
      }
      #Get UTR3hash of the current mRNA studied if exist
      if (exists $UTR3->{$mRNAname}){
        my %UTR3Hash=@{$UTR3->{$mRNAname}};
        brickToprintHash($output, \%UTR3Hash);
      }
    }
  } # tRNA
  elsif(exists  $tRNA->{$geneName}){
    my %tRNAHash=@{$tRNA->{$geneName}};
    foreach my $tRNAnum (keys %tRNAHash) {
      my $current_tRNAHash=$tRNAHash{$tRNAnum};
      brickToprintTab($output, $current_tRNAHash);
      my $tRNAname=@{$current_tRNAHash}[0];

      my %exonHash=@{$exon->{$tRNAname}};   
      brickToprintHash($output, \%exonHash);
      #Get cdshash of the current tRNA studied
      if (exists $cds->{$tRNAname}){
        my %cdsHash=@{$cds->{$tRNAname}};
        brickToprintHash($output, \%cdsHash);
      }
      #Get UTR5hash of the current tRNA studied if exist
      if (exists $UTR5->{$tRNAname}){
        my %UTR5Hash=@{$UTR5->{$tRNAname}};
        brickToprintHash($output, \%UTR5Hash);
      }
      #Get UTR3hash of the current tRNA studied if exist
      if (exists $UTR3->{$tRNAname}){
        my %UTR3Hash=@{$UTR3->{$tRNAname}};
        brickToprintHash($output, \%UTR3Hash);
      }
    }
  }   
}

#This function sort hash to give key in order from one to ...
sub sortByPos {
  my ($hashRef) = @_;
  foreach my $key (keys %{$hashRef}) { 
    my %hashDeep = @{ ${$hashRef} {$key} };

    my $cpt=1;
    my %hashtmp;
      foreach my $kDeep ( sort ({ $hashDeep{$a}[4] <=> $hashDeep{$b}[4] } keys %hashDeep)){  
        @{ $hashtmp{$cpt} } = @{$hashDeep{$kDeep}};
        $cpt++;       
    }
    @{ ${$hashRef} {$key} } =  %hashtmp ;
  }
} 

sub parseGFF {
  
  my($file_in,$fileName, $opt_expand) = @_;
  print( "Reading features from $fileName...\n");

  my %genes; my %trnagenes; my %mRNA;my %exon; my %cds; my %UTR5; my %UTR3; my %tRNA; my %nctpss; my %ncfpss; my %pieceStudied;
  my $cptgenes; my $cptmRNA;my $cptexon; my $cptcds; my $cptUTR5; my $cptUTR3; my $cpttRNA; my $cptgenetRNA; my $cptnctpss; my $cptncfpss;

  # read file and decompose it
  while (my $feature = $file_in->next_feature() ) {
      my $seqname = $feature->seq_id();#print "$seqname \n";
      my $source = $feature->source_tag();#print "source =$source \n";
      my $type = $feature->primary_tag();#print "type= $type \n";
     #Manage feature position // shoud be always sorted
      my $start = $feature->start();
      my $end = $feature->end();
      if ($start > $end){my $tmp=$start; $start=$end; $end=$tmp;} 
      my $score =  $feature->score();if (! defined $score){$score = ".";}
      my $strand  = $feature->strand();if ($strand eq "1"){$strand="+";}elsif($strand eq "0"){$strand="-";}#print "strand= $strand \n";
      my $frame  = $feature->frame(); if ( ! defined $frame){$frame = ".";} #print "frame= $frame";
      
      my @ID = $feature->get_tag_values('ID');
      my $ID=$ID[0];#print "ID=  $ID\n";
      my $ParentID;

      my @groups;
      my @tags = $feature->get_all_tags(); 
      foreach my $tag  (@tags){
        my @tagValue = $feature->get_tag_values($tag);
        push (@groups, [$tag, $tagValue[0]]);
        if ($tag eq "Parent"){$ParentID=$tagValue[0];}
      }

      if ( $type eq 'gene' ) {
        if ($ID =~ m/trnascan/){
          $cptgenetRNA++;
        push( @{ $trnagenes{$ID}}, ($ID,$seqname,$source,$type,$start,$end,$score,$strand,$frame,\@groups) );
        }
        else{   
        my $pieceStudiedName="$seqname$strand";
        push( @{ $pieceStudied{$pieceStudiedName}}, $ID );
        push( @{ $genes{$ID}}, ($ID,$seqname,$source,$type,$start,$end,$score,$strand,$frame,\@groups) );
        $cptgenes++;
        }
      }
      else {
#     $output->write_feature($feature);     
        if ($type eq 'mRNA'){
          my %mRNAsInfo; $cptmRNA++;
          if ( !defined( $genes{$ParentID}) ) {
            printf( STDERR "mRNA -> gene parent ".$ParentID." dont exists. Line= $ID, $seqname,$source,$type,$start,$end,$score,$strand,$frame,\@groups\n"); exit();
          }
          if (exists $mRNA{$ParentID}){
              # get the hash      
              my %mRNAsInfo=@{ $mRNA{$ParentID}};
              #nb element in hash
              my $nbr = keys (%mRNAsInfo);
       #       printf ("\nthere is $nbr keys described\n");
              # add value in hash
              push ( @{ $mRNAsInfo{$nbr+1} }, ($ID,$seqname,$source,$type,$start,$end,$score,$strand,$frame,\@groups) );
              # put hash in hash
              @{ $mRNA{$ParentID} } =  %mRNAsInfo ;    
            }
            else{
              @{ $mRNAsInfo{1} }= ($ID,$seqname,$source,$type,$start,$end,$score,$strand,$frame,\@groups);
              push ( @{ $mRNA{$ParentID} },  %mRNAsInfo);

            }
          
        }

        elsif ($type eq "exon"){
          $cptexon++;
#<<<<<<<<<<<<<<<< Expand exon If seceral parents >>>>>>>>>>>>>>>
            #We will expand exon if necessary
            my @Parent = $feature->get_tag_values('Parent');
            my $sizeParent=$#Parent;
            # for each parent of the exon, create an exon (even if only one...)
            my $nbmRNAexpanded=0;
            foreach my $mRNA_ID (@Parent){
              $nbmRNAexpanded++;       
              my @groupsR;
              if ($sizeParent != 0){ # If this exon need to be expand
                      if ($nbmRNAexpanded != 1){ #dont count the original exon
                        $nbexpand++;
                      }
                     # print "expand mRNA $mRNA_ID\n";
                      @tags = $feature->get_all_tags(); 
                      foreach my $tag  (@tags){
                        if ($tag eq "Parent"){
                          push (@groupsR, [$tag, $mRNA_ID]);
  #                        print " test1=$mRNA_ID\n";
                        }
                       elsif($tag eq "ID"){ # Change name of expanded exon tin order to have a unique identifier/ID
                          my @tagValue = $feature->get_tag_values($tag);
                          push (@groupsR, [$tag, "$tagValue[0]]-$mRNA_ID"]);
                       }
                       else{
                         my @tagValue = $feature->get_tag_values($tag);
                         push (@groupsR, [$tag, $tagValue[0]]);
  #                       print "tagvalue = $tagValue[0] "
                        }
                       
                      }
              }
              else {@groupsR=@groups;} # If this exon cannot be expand

              #Now save result in Hash
              my %exonsInfo=();
              #Test if exon hash already exists
              if (exists $exon{$mRNA_ID}){
                # get the hash      
                my %exonsInfo=@{ $exon{$mRNA_ID}};
                #nb element in hash
                my $nbr = keys (%exonsInfo);
                # add value in hash
                push ( @{ $exonsInfo{$nbr+1} }, ($ID,$seqname,$source,$type,$start,$end,$score,$strand,$frame,\@groupsR) );
                # put hash in hash
                @{ $exon{$mRNA_ID} } =  %exonsInfo ;
              }
              else{
                @{ $exonsInfo{1} }= ($ID,$seqname,$source,$type,$start,$end,$score,$strand,$frame,\@groupsR);           
                push ( @{ $exon{$mRNA_ID} },  %exonsInfo);
              }
            } # end foreach (even if only one loop)
        }

        elsif ($type eq "CDS"){
          $cptcds++;
          my %cdsInfo;
          my $mRNA_ID=$ParentID;
          #Test if cds hash already exists
          if (exists $cds{$mRNA_ID}){
            # get the hash      
            my %cdsInfo=@{ $cds{$mRNA_ID}};
            #nb element in hash
            my $nbr = keys (%cdsInfo);
    #        printf ("\nthere is $nbr keys described\n");
            # add value in hash
            push ( @{ $cdsInfo{$nbr+1} }, ($ID,$seqname,$source,$type,$start,$end,$score,$strand,$frame,\@groups) );
            # put hash in hash
            @{ $cds{$mRNA_ID} } =  %cdsInfo ;    
          }
          else{
            @{ $cdsInfo{1} }= ($ID,$seqname,$source,$type,$start,$end,$score,$strand,$frame,\@groups);
            push ( @{ $cds{$mRNA_ID} },  %cdsInfo);
          }

        }

        elsif ($type eq "five_prime_UTR"){
          $cptUTR5++;
          my %UTR5Info;
          my $mRNA_ID=$ParentID;
          if (exists $UTR5{$mRNA_ID}){
            # get the hash      
            my %UTR5Info=@{ $UTR5{$mRNA_ID}};
            #nb element in hash
            my $nbr = keys (%UTR5Info);
            # add value in hash
            push ( @{ $UTR5Info{$nbr+1} }, ($ID,$seqname,$source,$type,$start,$end,$score,$strand,$frame,\@groups) );
            # put hash in hash
            @{ $UTR5{$mRNA_ID} } =  %UTR5Info ;    
          }
          else{
            @{ $UTR5Info{1} }= ($ID,$seqname,$source,$type,$start,$end,$score,$strand,$frame,\@groups);
            push ( @{ $UTR5{$mRNA_ID} },  %UTR5Info);
          }
        }

        elsif ($type eq "three_prime_UTR"){
          $cptUTR3++;
          my %UTR3Info;
          my $mRNA_ID=$ParentID;
          if (exists $UTR3{$mRNA_ID}){
            # get the hash      
            my %UTR3Info=@{ $UTR3{$mRNA_ID}};
            #nb element in hash
            my $nbr = keys (%UTR3Info);
            # add value in hash
            push ( @{ $UTR3Info{$nbr+1} }, ($ID,$seqname,$source,$type,$start,$end,$score,$strand,$frame,\@groups) );
            # put hash in hash
            @{ $UTR3{$mRNA_ID} } =  %UTR3Info ;    
          }
          else{
            @{ $UTR3Info{1} }= ($ID,$seqname,$source,$type,$start,$end,$score,$strand,$frame,\@groups);
            push ( @{ $UTR3{$mRNA_ID} },  %UTR3Info);
          }
        }

        elsif($type eq "tRNA"){
          $cpttRNA++;
          my %tRNAsInfo;
          my $geneID = $ParentID;
          if ( !defined( $trnagenes{$geneID}) ) {
            printf( STDERR "tRNA -> geneID ".$geneID." dont exists. Line= $ID, $seqname,$source,$type,$start,$end,$score,$strand,$frame,\@groups\n"); exit();
          }
          if (exists $tRNA{$geneID}){
              # get the hash      
              my %tRNAsInfo=@{$tRNA{$geneID}};
              #nb element in hash
              my $nbr = keys (%tRNAsInfo);
              # add value in hash
              push ( @{ $tRNAsInfo{$nbr+1} }, ($ID,$seqname,$source,$type,$start,$end,$score,$strand,$frame,\@groups) );
              # put hash in hash
              @{ $tRNA{$geneID} } =  %tRNAsInfo ;    
            }
            else{
              @{ $tRNAsInfo{1} }= ($ID,$seqname,$source,$type,$start,$end,$score,$strand,$frame,\@groups);
              push ( @{ $tRNA{$geneID} },  %tRNAsInfo);

            }
        }
        elsif($type eq "non_canonical_three_prime_splice_site"){
          my %spliceInfo;$cptnctpss++;
          if (exists $nctpss{$ParentID}){
              # get the hash      
              my %spliceInfo=@{$nctpss{$ParentID}};
              #nb element in hash
              my $nbr = keys (%spliceInfo);
              # add value in hash
              push ( @{ $spliceInfo{$nbr+1} }, ($ID,$seqname,$source,$type,$start,$end,$score,$strand,$frame,\@groups) );
              # put hash in hash
              @{ $nctpss{$ParentID} } =  %spliceInfo ;    
            }
            else{
              @{ $spliceInfo{1} }= ($ID,$seqname,$source,$type,$start,$end,$score,$strand,$frame,\@groups);
              push ( @{ $nctpss{$ParentID} },  %spliceInfo);

            }
        }
        elsif($type =~ m/non_canonical_five_/){
          my %spliceInfo;$cptncfpss++;
          if (exists $ncfpss{$ParentID}){
              # get the hash      
              my %spliceInfo=@{$ncfpss{$ParentID}};
              #nb element in hash
              my $nbr = keys (%spliceInfo);
              # add value in hash
              push ( @{ $spliceInfo{$nbr+1} }, ($ID,$seqname,$source,$type,$start,$end,$score,$strand,$frame,\@groups) );
              # put hash in hash
              @{ $ncfpss{$ParentID} } =  %spliceInfo ;    
            }
            else{
              @{ $spliceInfo{1} }= ($ID,$seqname,$source,$type,$start,$end,$score,$strand,$frame,\@groups);
              push ( @{ $ncfpss{$ParentID} },  %spliceInfo);

            }
        }

        else{
            printf( STDERR "Feature skipped: $type... Use option\n");
            exit();
        }
      }
  }
  
  $ref_istream->close();
  printf( "Read %d genes.(tRNA included if present)\n", $cptgenes );
  my $nbGeneFrommRNA = keys (%mRNA);
  printf( STDERR "Read genes $nbGeneFrommRNA (tRNA excluded if present)\n");
  my $nbmRNA = keys (%exon);
  printf( STDERR "Read mRNAs $nbmRNA\n");

#  if ($opt_expand eq "yes"){ print "$nbexpand exons added during the expansion\n";}

  return \%genes, \%trnagenes, \%mRNA, \%exon, \%cds, \%UTR5, \%UTR3, \%tRNA, \%pieceStudied;
}


__END__

=head1 NAME

maker_checkFusionSplitBetweenTwoBuilds.pl - Compare two gene build in GFF3 format in order 
to detect the gene fron build 1 (--ref file) that are split or fused in the gene build 2 (--tar file). The result is written in the file "outputSplitMergeRef" (For gene involved from build1) and
"outputSplitMergeTar" (For gene involved from build2) in the specified output directory. In more, complex cases 
where more than 2 genes of each build are overlapping, we sort them in a separate file call Cluster as well written in the specidied output directory. 

=head1 SYNOPSIS

    ./maker_checkFusionSplitBetweenTwoBuilds.pl --ref=infile --tar=infile --output=outDirectory
    ./maker_checkFusionSplitBetweenTwoBuilds.pl --help

=head1 OPTIONS

=over 8

=item B<--ref>, B<--reffile> or B<-f>

Input GFF3 file correponding to gene build 1.

=item B<--tar>, B<--tarfile> or B<-t>

Input GFF3 file corresponding to gene build 2.

=item  B<--out>, B<--output> or B<-o>

Output directory where diffrent output files will be written.

=item B<--help> or B<-h>

Display this helpful text.

=back

=cut
