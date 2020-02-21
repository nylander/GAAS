#!/usr/bin/env perl

##########################
# Jacques Dainat 11/2015 # 
###########################

#libraries
use File::Basename;
use strict;
use warnings;
use Data::Dumper;
use Carp;
use Bio::SimpleAlign;
use Bio::LocatableSeq;
use Getopt::Long;
use IO::File;
use Pod::Usage;
use Bio::SeqIO;
use List::MoreUtils qw(uniq);
use Storable 'dclone';
use Clone 'clone';

# END libraries
# PARAMETERS - OPTION
my $opt_position;
my $opt_output=undef;
my $opt_nbByWindow=undef;
my $opt_windowsSize=undef;
my $opt_aliDir=undef;
my $opt_tsv=undef;
my $opt_cons=undef;
my $opt_consList=undef;
my $opt_consThreshold=undef;

my $opt_help = 0;
# END PARAMETERS - OPTION


# OPTION MANAGMENT
if ( !GetOptions( 'position|p=s' => \$opt_position,
                  'tsv=s'      => \$opt_tsv,
                  'o|output=s'      => \$opt_output,
                  'v|value|window=i'      => \$opt_windowsSize,
                  'ad=s'      => \$opt_aliDir,
                  'consensus|c' => \$opt_cons,
                  'consensus_list|cl=s' => \$opt_consList,
                  'consensus_threshold|ct=i' => \$opt_consThreshold,
                  'nbbw=i'      => \$opt_nbByWindow,
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

if ((!$opt_position or !$opt_aliDir or !$opt_tsv)){ 
    pod2usage( {
           -message => "\nIf you want to merge sequence by window, at least 3 parameter is mandatory: position file, tsv file, and directory containing alignments\n".
          "Many optional parameters are available. Look at the help documentation to know more.\n",
           -verbose => 0,
           -exitval => 2 } );
}
#############################
####### Manage options #######
#############################

#Choose y axis value:
if(! $opt_nbByWindow){
  $opt_nbByWindow="4";
  print "you didn't choose any value for the number of minimum gene that must be present to keep a window. We will take the default value.\n";
}
print "Number minimum of gene by window : $opt_nbByWindow .\n";
#Choose opt_windowsSize value:
if(! $opt_windowsSize){
  print "you didn't choose any value for the window. We will take the default value.\n";
  $opt_windowsSize=1;
}
print "Window size used : $opt_windowsSize Mb.\n";

if($opt_cons ){
  print "You decided to create consensus sequences in the alignements when several individual from the same species are present.\n";
}

my %consensListOk;
if($opt_consList){
  print "Moreover You decided to make consensus between species. Here is the detail:\n";
  # Mange consensus list given 
  my @consensListPair;

  @consensListPair= split(/,/, $opt_consList);
  my $nbCons=$#consensListPair+1;
  print "You defined $nbCons species for creating consensus sequences.\n";
  foreach my $consensusTuple (@consensListPair){
    my @consensList= split(/\//, $consensusTuple);
    if($#consensList != 1){ # Attribute alone
      print "Problem with this tuple $consensusTuple. You should have something like: Species/consensusName\n";
    }
    else{ # Attribute we have to replace by a new name
      push(@{$consensListOk{$consensList[1]}}, $consensList[0]);
    }
  }
  foreach my $key (keys %consensListOk){ 
    print "You dicided to make a consensus of (inter)species <@{$consensListOk{$key}}> and call this new sequence <$key>\n";
  }
}

if($opt_cons or $opt_consList){
  if(! $opt_consThreshold){
    $opt_consThreshold=0;
  }
  print "The consensus residue has to appear at least threshold $opt_consThreshold% of the sequences at a given location, otherwise a '?'\n";
}

#############################
####### Manage output #######
#############################

my $report_output;

if( $opt_output){

  $opt_output=~ s/\..*//g;
  $report_output=$opt_output."_report";
}
else{
  $opt_output="output_prepare_matrice_by_window";
  $report_output=$opt_output."_report";;
}

if (-d $opt_output ){
    print "Directory $opt_output already exists.\n";exit;
}
else{
    mkdir $opt_output;
    mkdir $report_output;
    print "$opt_output dir created\n";
}

                            #######################
                            #        MAIN         #
                            #######################

########## constant #############

my $nbFile=0;
my @listTmpFile;

#######################
# Manage input files  #
#######################

### READ POSITION FILE AND STORE RESULTS
my $position_stream = IO::File->new();
$position_stream->open( $opt_position, 'r' ) or
croak(
     sprintf( "Can not open '%s' for reading: %s", $opt_position, $! ) );

my %hash_chr2name;
my $count=0;
  foreach my $line (<$position_stream>){
    chomp($line);
    if($count >=1){ # avoid header
      my @list=split /@/,$line;
      #geneName@Chromosome@position@GC

      $hash_chr2name{$list[1]}{$list[2]}= $list[0];
    }
    $count++;
  }

### READ TSV FILE AND STORE RESULTS
my $tsv_stream = IO::File->new();
$tsv_stream->open( $opt_tsv, 'r' ) or
croak(
     sprintf( "Can not open '%s' for reading: %s", $opt_tsv, $! ) );

my %hash_HordeumName2Aliname;

  foreach my $line (<$tsv_stream>){
    chomp($line);
    my @list=split /\t/, $line;
    #print @list."\n";
    #MarqueurName Nbspecies NbIndividu  NbSequence  SequenceNamesList(separated by tabulation)
    my $count=0;
    my $marqueurName;
    my $hordeumName;
    foreach my $el (@list){

      if($count == 0){ # marqueurName
        $marqueurName=$el;
      }
      if($count >=4){ # SequenceNamesList

        if ($el =~ "H_vulgare"){ # H_vulgare sequence => shoud contain EPlHVUG or MLOC that will be used to retrieve position information
         
          if($el =~ "EPlHVUG"){
            
            $el =~ /.*_(EPlHVUG[^|]*)\|.*/;
            $hordeumName = $1;
          }
          else{
            $el =~ /.*\|([^\.]*).*/;
            $hordeumName = $1;
            last;
          }
          
        }
      }
    $count++;
    }

    $hash_HordeumName2Aliname{$hordeumName}=$marqueurName;
  }


### READ ALI FOLDER ####
my @listFile;
if( -d $opt_aliDir){
  opendir(DIR, $opt_aliDir) or die "cannot open dir $opt_aliDir: $!";
  @listFile= readdir DIR;
  closedir DIR;
}
else{
  print "$opt_aliDir is not a directory or doesn't exists !\n";exit;
}

##########################
#    General results     #
##########################


my $currentCenter;
my $path;
my %Ali;my $tmpAli=\%Ali;


#For each chromosome
foreach my $chr (keys %hash_chr2name){

  print "Process chromosome".$chr."\n";
  my $nbWinOk=0;
  my $dir=$opt_output."/".$chr;
  mkdir $dir;
  my $currentLimit = 0;

  #report
  my $report_stream = IO::File->new();
  $report_stream->open( $report_output."/".$chr, 'w' ) or croak( sprintf( "Can not open '%s' for reading: %s", $report_output, $! ) );

  # need to know the highest value
  my $highvalue=get_highest_value(\%hash_chr2name, $chr);

  
  # Go through all the chromosome window by window
  my $stop=0;
  while ($stop != 2){ 

    if( (1000000*$currentLimit) > $highvalue){$stop++};
    
    my $nbPrintedVal=0;
    $currentLimit = $currentLimit + $opt_windowsSize;
    $currentCenter=$currentLimit-($opt_windowsSize/2);
    $currentCenter=sprintf "%.1f",$currentCenter;
    $path=$dir."/".$currentCenter."Mb.aln";
    #print "study of chr $chr window ".1000000*($currentLimit-$opt_windowsSize)."-".1000000*$currentLimit."\n";

    %$tmpAli = (); # empty the hash

    # we go through the chromosome positions by increasing order
    foreach my $position (sort {$a<=>$b} keys $hash_chr2name{$chr}){
      #print "study of $hash_chr2name{$chr}{$position} $position\n";
      if($position > (1000000*($currentLimit-$opt_windowsSize)) and $position < (1000000*$currentLimit) ){ #kb to nucleotide : * 1000000

        if(exists ($hash_HordeumName2Aliname{$hash_chr2name{$chr}{$position}})){
          
         if (fill_for_ali($tmpAli, $hash_HordeumName2Aliname{$hash_chr2name{$chr}{$position}}, $nbPrintedVal)){ # save sequence of the alignement to create supermatrice
          $nbPrintedVal++;
         }
        }
        else{
          #print "Hordeum Name: $hash_chr2name{$chr}{$position} doesn't exists among the Hordeum names of alignments we are working with.\n";
        }
          
      }
      # We are over the chromosome size, consequently we stop
      elsif($position > (1000000*($currentLimit-$opt_windowsSize))){last;}
    }
    
    #report
    print $report_stream $currentCenter."Mb\t".$nbPrintedVal."\n";

    # remove files if not enough value  ín the interval
    if($nbPrintedVal < $opt_nbByWindow){
      #print "only $nbPrintedVal alignment ! We delete it \n";
      unlink $path;
    }
    else{ # enough value in the window, so we print the result
      $nbWinOk++;
      my $ostreamChrFile  =  Bio::SeqIO->new(-format => "fasta", -file => ">$path");

      ###################
      #Create consensus
      my %listOfAli;
      if($opt_cons){
        my %hashNum;
        my %hashNames;
        #create list name of the current Ali
        foreach my $key (keys %$tmpAli ){
          #print "name = $key\n";
          $key =~ /(.*_[^_]*)_/;
          my $nameRough = $1;
          $hashNum{$nameRough}++;
          push(@{$hashNames{$nameRough}}, $key); 
        }

        #
        foreach my $nameRough (keys %hashNum ){
          if($hashNum{$nameRough} > 1){          
            # create a small alignementof only same species
            my $aln = new Bio::SimpleAlign();     
            foreach my $nameComplete (@{$hashNames{$nameRough}}){
              my $seq = new Bio::LocatableSeq(-seq => $tmpAli->{$nameComplete} , -id  => "$nameComplete");
              $aln->add_seq($seq);
              print "We put $nameComplete in $nameRough \n";
              #remove the sequence from old ali no more usefull
              delete $tmpAli->{$nameComplete}; #REMOVE THAT KEY-VALUE pair
            }
            #print "we keep name $nameRough -------\n";
            $listOfAli{$nameRough}=$aln;
          }
        }
      }

      if($opt_consList){
        foreach my $ConsensName (keys %consensListOk ){
          
          my $SuperAln = undef;
          foreach my $nameCons (@{$consensListOk{$ConsensName}} ){
            #print "Does it start by ? $nameCons\n";
            #foreach my $na ( keys %$tmpAli){
            #  print "$na \n";
            #}
            # check in the original alignement
            if (grep {/$nameCons.*/} keys %$tmpAli){

              if(! $SuperAln){ $SuperAln = new Bio::SimpleAlign();}

              foreach my $name (keys %$tmpAli){
                if($name =~ /$nameCons.*/){
                  my $seq = new Bio::LocatableSeq(-seq => $tmpAli->{$name} , -id  => "$name");
                  $SuperAln->add_seq($seq);
                  print "We put $name in $ConsensName \n";
                  #remove the sequence from old ali no more usefull
                  delete $tmpAli->{$name}; #REMOVE THAT KEY-VALUE pair
                }
              }
            }

            #check if alignement that can have been done by the option "c"
            if (grep {/$nameCons.*/} keys %listOfAli){

              if(! $SuperAln){ $SuperAln = new Bio::SimpleAlign();}

              my %deep_copy_listOfAli = %{ clone (\%listOfAli) };
              foreach my $key (%deep_copy_listOfAli){
                if ($key =~ /$nameCons.*/){

                  my $currentAli=$listOfAli{$key};

                  foreach my $seq ($currentAli->each_seq) {
                    $SuperAln->add_seq($seq);
                  }
                  print "We put $key in $ConsensName \n";
                  delete $listOfAli{$key}; #REMOVE THAT KEY-VALUE pair
                } 
              }
            }
          }
          if($SuperAln){
            $listOfAli{$ConsensName}=$SuperAln;
          }
        }
      }

      #create a consensus from each ali if asked
      foreach my $consensName (keys %listOfAli){ 
        my $aln=$listOfAli{$consensName};
        my $consensus = $aln->consensus_string($opt_consThreshold);
        $consensus=~ s/\?/-/g; #  quand il y a que des gap il remplace par ? .Donc ici on remplace ? par -

        #now add the consensus sequence to
        $tmpAli->{$consensName."_Cons"}=$consensus;
      }

      # End create consensus
      ######################

      #print result alignment
      foreach my $key (keys %$tmpAli ){

        my $seq = Bio::Seq->new(-display_id => $key, -seq => $tmpAli->{$key});
        $ostreamChrFile->write_seq($seq);
      }

      #print "$nbPrintedVal ali read\n";
      $position_stream->close();
    }
  }

  # remove folder in no file
  if (is_folder_empty($dir)) {
    rmdir $dir
  }
  else{
    print "We found $nbWinOk window(s) of $opt_windowsSize Mb with more than $opt_nbByWindow alignments\n";
  }
}
$position_stream->close();
print("Parsing Finished\n\n");


      ######################### 
      ######### END ###########
      #########################
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

sub fill_for_ali{
  my ($hash,$name,$nbPrintedVal)=@_;

  ### get file name !! 
  my @matches = grep { /$name/ } @listFile;
  if (@matches >= 2){
    print "several file start with the name $name whithin the directory $opt_aliDir. We don't know which one take. We stop. \n";
    return 0;
  }
  if (@matches == 0){
    #print "No file with name starting by $name found whithin the directory $opt_aliDir. We don't know which one take. We stop. \n";
    return 0;
  }
  else{
    $name=$matches[0];
    
    if(-f "$opt_aliDir/$name"){

      #print "$opt_aliDir/$name\n";
      my $seqio = Bio::SeqIO->new(-file => "$opt_aliDir/$name", '-format' => 'Fasta');
      my %idPerformed;
      my $string;

      #get actualSize
      my $sizeAlAli=0;
      foreach my $key (keys %$hash){
        $sizeAlAli = length($hash->{$key});
        last;
      }

      #print "read file $opt_aliDir/$name\n";
      while(my $seq = $seqio->next_seq) {

        $string = $seq->seq;
        my $id_original = $seq->primary_id;
        my @ids= split /_/,$id_original ;
        my $newID=$ids[0]."_".$ids[1]."_".$ids[2];
        #print "newID $newID\n";
        
        if(! exists($hash->{$newID})){ # If Id is a new one
          if(! $nbPrintedVal == 0){ # if not the first alignment read we have to add empty seq in front 
            
            # If it's new but some other were already existing we add empty seq in front
            my $emptySeq = "";
            $emptySeq =~ s/^(.*)/'-' x $sizeAlAli . $1/mge; # create a string of gap for the sequence to add in front of the ali
            my $NewSeq="$emptySeq$string" ;
            $hash->{$newID}=$NewSeq;
          }
          else{ #first ali file
            #print "save $newID\n";
            $hash->{$newID}=$string;
          }
        }
        else{ # If Id already known
          if(! exists ($idPerformed{$newID})){ #If we have not already seen this ID this round ! Allow avoiding duplication !!
            #print "already exists we append it $newID\n";
            $hash->{$newID}.=$string;
          }
          else{print "Duplication removed !\n";}
        }
        $idPerformed{$newID}++;
      }
      # check if we have to add empty seq (sequences present in pevious alignment that are not present in the current one)
      foreach my $key (keys %$hash){
        if (! exists ($idPerformed{$key})){ # add empty line for that seuquence that doesent exists
          my $n = length($string);
          my $emptySeq = "";
          $emptySeq =~ s/^(.*)/'-' x $n . $1/mge;
          $hash->{$key}.=$emptySeq;
        }
      }

      return 1;
    }
    else{print "File $name not present within the folder $opt_aliDir...\n" ;
         return 0;}
  }
}

sub get_highest_value{

  my ($hash,$chr) = @_;
  my $high=0;

  foreach my $position (keys $hash->{$chr}){
    if($position > $high){
      $high=$position ;
    }
  }
  return $high;
}

sub is_folder_empty {
    my $dirname = shift;
    opendir(my $dh, $dirname) or die "Not a directory";
    return scalar(grep { $_ ne "." && $_ ne ".." } readdir($dh)) == 0;
}

__END__

=head1 NAME

script.pl -
The script take a position file as input and a directory containing the alignments. It will concatenates files to create mini-superMatrices corresponding to the size of window size choosen.
The aim is to check if ther is a clear phylogeny pattern by portion of the chromosome. If hte pattern is different in different location of hte chromosome it could be a hint of recombination.

=head1 SYNOPSIS

    ./script.pl -p positionFile --ad directoryWithAlignments [ --output outfile ]
    ./script.pl --help

    ./prepare_matrice_by_window.pl -p Hordeum_position.txt -tsv clusters_1L_7speIN_1speOUT_noDup_withHordeum6235.tsv -ad ALIGNMENTS -c -o out_consens

    example with -cl option
    ./prepare_matrice_by_window.pl -p Hordeum_position.txt -tsv clusters_1L_7speIN_1speOUT_noDup_withHordeum6235.tsv -ad aliTest/ -c -cl Ae_longissima/D,Ae_sharonensis/D,Ae_bicornis/D,Ae_searsii/D,Ae_tauschii/D,Ae_comosa/D,Ae_uniaristata/D,Ae_umbellulata/D,Ae_caudata/D,T_urartu/A,T_boeoticum/A,Ae_speltoides/B,Ae_mutica/B

=head1 OPTIONS

=over 8

=item B<-p> or B<--position>

Input file containing position information related to a reference in that format:
GeneName@ChromosomeName@PositionInNucleotide

=item B<--tsv>
This option define the tsv file that contain information needed to link the position to the alignments. Format like this:
Nom marqueur (sequence untiliser pour construire le cluster) \t nombre sp \t nombre individu \t nombre sequences \t Toutes les entetes des sequences de l'alignement seaparer par une tabulation 

=item B<--ad>

Input directory containing alignment files in fasta format.  Name of files should be the same those find at "GeneName" from the position file.

=item B<-v>, B<--value> or B<--window>

Allows to define the sliding window size to use to clusterize the alignments (in megabase) (Default 1).

=item B<--nbbw>

Number minimum of sequence that must be present within a window to take it into account (Default 4).

=item B<--consensus>, B<-c>

This option when activated, will create consensus sequences in the alignements when several individual from the same species are present.
Ae_comosa_Tr272 => Using the two first word linked by an underscore to define the species name. Here Ae_comosa will be used.

=item B<--consensus_list>, B<-cl>
This option when activated, will create consensus sequences in the alignements of sequences starting by names defined, and will be grouped by Names chosen. Explantion below:
Example: Species1/A,Species2/A,Species3/B,Species4/B
In this example Species1 and Species2 will create a consensus called A, and species3 and species4 will create a consensus called B.

=item B<--consensus_threshold>, B<-ct>
Optional treshold ranging from 0 to 100. The consensus residue has to appear at least threshold % of the sequences at a given location, otherwise a '?' character will be placed at that location.
Default value=0;

=item B<-o> or B<--output> 

Output name of the directory that will contain results

=item B<-h> or B<--help>

Display this helpful text.

=back

=cut