#!/usr/bin/env perl

use strict;
use warnings;
use POSIX qw(strftime);
use File::Basename;
use Carp;
use Getopt::Long;
use IO::File;
use Pod::Usage;
use Statistics::R;
use BILS::CheckModule qw(:Ok);
use BILS::Plot::R qw(:Ok);
use BILS::Handler::GFF3handler qw(:Ok);
use BILS::Handler::GXFhandler qw(:Ok);
use Bio::Tools::GFF;
use BILS::GFF3::Statistics qw(:Ok);

my $header = qq{
########################################################
# BILS 2016 - Sweden                                   #  
# jacques.dainat\@bils.se                               #
# Please cite BILS (www.bils.se) when using this tool. #
########################################################
};


my @opt_files;
my $opt_output=undef;
my $opt_breaks;
my $Xpercent=1;
my $opt_help = 0;

my @copyARGV=@ARGV;
if ( !GetOptions( 'f|gff|ref|reffile=s' => \@opt_files,
                  'o|out|output=s' => \$opt_output,
                  'w|window|b|break|breaks=i'      => \$opt_breaks,
                  'x|p=i'      => \$Xpercent,
                  'h|help!'         => \$opt_help ) )
{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

if ($opt_help) {
    pod2usage( { -verbose => 2,
                 -exitval => 2,
                 -message => "$header\n" } );
}

if ( ! defined( $#opt_files  >= 0) ) {
    pod2usage( {
           -message => "$header\nMust specify at least 1 parameters:\nReference data gff3 file (--gff)\n",
           -verbose => 0,
           -exitval => 1 } );
}

# #######################
# # START Manage Option #
# #######################

if (defined($opt_output) ) {
  if (-d $opt_output){
    print "The output directory choosen already exists. Please geve me another Name.\n";exit();
  }
  else{
    mkdir $opt_output;
  }
}

my $ostreamReport;
if (defined($opt_output) ) {
  $ostreamReport=IO::File->new(">".$opt_output."/report.txt" ) or croak( sprintf( "Can not open '%s' for writing %s", $opt_output."/report.txt", $! ));
}
else{
  $ostreamReport = \*STDOUT or die ( sprintf( "Can not open '%s' for writing %s", "STDOUT", $! ));
}
my $string1 = strftime "%m/%d/%Y at %Hh%Mm%Ss", localtime;
$string1 .= "\n\nusage: $0 @copyARGV\n\n";

print $ostreamReport $string1;
if($opt_output){print $string1;}


#############################
####### Manage R option #####
#############################

#Choose breaks value:
if(! $opt_breaks){
  $opt_breaks="1000";
}

#############################
####### Manage output #######
#############################
my $outputPDF;
if (defined($opt_output) ) {
  if (-f $opt_output){
      print "Cannot create a directory with the name $opt_output because a file with this name already exists.\n";exit();
  }
 $outputPDF=$opt_output."/intronPlot.pdf";   
}
else{
  $outputPDF="intronPlot.pdf";
}

# Check R is available. If not we try to load it through Module software
if ( system("R --version 1>/dev/null 2>/dev/null") == 0 ) {
  print "R is available. We can continue\n";
}
else {
  print "R is not loaded. We try to load it.\n";
  if(module_software_installed){
    module_load("R");
  }
  else{
    print "Module tool doesn't exists. We cannot load R through it.";
  }
}

# #####################################
# # END Manage OPTION  
# #####################################



                                                      #######################
                                                      #        MAIN         #
#                     >>>>>>>>>>>>>>>>>>>>>>>>>       #######################       <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#PART 1
###################################
# Read input gff3 files one by one and save value in hash of list

my @introns;
foreach my $file (@opt_files){

  print "Reading ".$file,"\n";

  # parse file name te remove extension
  my ($file1,$dir1,$ext1) = fileparse($file, qr/\.[^.]*/);

  ######################
  ### Parse GFF input #
  my ($hash_omniscient, $hash_mRNAGeneLink) = BILS::Handler::GFF3handler->slurp_gff3_file_JD($file);
  print("Parsing Finished\n\n");
  ### END Parse GFF input #
  #########################

  #print statistics
  my $stat = gff3_statistics($hash_omniscient);
  #print statistics
  foreach my $infoList (@$stat){
    foreach my $info (@$infoList){
      print $ostreamReport "$info";
    }
    print $ostreamReport "\n";
  }


  ######################
  ### Parse GFF input #
  # get nb of each feature in omniscient;
  foreach my $tag_l2 (keys %{$hash_omniscient->{'level2'}}){
    foreach my $id_l1 (keys %{$hash_omniscient->{'level2'}{$tag_l2}}){
      my $one_f2 = $hash_omniscient->{'level2'}{$tag_l2}{$id_l1}[0];

      #######################
      #get feature1 and info
      my $feature_l1=undef;
      my $tag_l1;
      foreach my $tag_level1 (keys %{$hash_omniscient->{'level1'}}){
        if (exists ($hash_omniscient->{'level1'}{$tag_level1}{$id_l1})){
          $feature_l1=$hash_omniscient->{'level1'}{$tag_level1}{$id_l1};
          $tag_l1=$tag_level1;
          last;
        }
      }
      if(! $feature_l1){print "Problem ! We didnt retrieve the level1 feature with id $id_l1\n";exit;}

      #####
      # get all level2
      my $All_l2_single=1;
      my $counterL2_match=-1;
      foreach my $feature_l2 ( @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}} ){

        #MATCH CASE - We ahve to count the L2 match features
        if($tag_l2 =~ "match"){
          my @sortedList = sort {$a->start <=> $b->start} @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}};
          my $indexLastL2 = $#{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}};
          $counterL2_match++;

          if($counterL2_match > 0 and $counterL2_match <= $indexLastL2){
            my $intronSize = $sortedList[$counterL2_match]->start - $sortedList[$counterL2_match-1]->end;
            push @introns, $intronSize;
          }
        }

        ######
        #get all level3
        my $id_l2=lc($feature_l2->_tag_value('ID'));

        foreach my $tag_l3 ( keys %{$hash_omniscient->{'level3'}} ){

          if(exists_keys($hash_omniscient,('level3',$tag_l3, $id_l2))){

          my $counterL3=-1;
          #Initialize intron to 0 to avoid error during printing results
          my $indexLast = $#{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2}};
          
          my @sortedList = sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2}};
          
            foreach my $feature_l3 ( @sortedList ){

              #count number feature of tag_l3 type
              $counterL3++;

              ################
              #Manage Introns# 
              # from the second intron to the last (from index 1 to last index of the table sortedList) 
              # We go inside this loop only if we have more than 1 feature.
              if($counterL3 > 0 and $counterL3 <= $indexLast){
                my $intronSize = $sortedList[$counterL3]->start - $sortedList[$counterL3-1]->end;
                push @introns, $intronSize;
              }
            }# END FOREACH L3
          }
        }
      }
    }
  }
}

# PART 2
###############################
my $biggest_value=0;
my $pathIntron="tmp_intron.txt";
my @sorted_intron = (sort { $a <=> $b } @introns);
#########################
# Write value in tmp files
 
# Manage Output
my $ostreamAED   = IO::File->new();
$ostreamAED->open( $pathIntron, 'w' ) or
      croak(
        sprintf( "Can not open '%s' for writing %s", $pathIntron, $! )
      );
foreach  my $value (@sorted_intron){
  print $ostreamAED "$value\n";
  if($value > $biggest_value){
    $biggest_value=$value;
  }
}
$ostreamAED->close();

# Part 3
#########################################
#Calcul longest after remove X percent  #
my $lastIndex = $#introns;
my $nbValueToRemove = int(($Xpercent*($lastIndex+1))/100);
my $resu =  $sorted_intron[$lastIndex-$nbValueToRemove];

my $stringPrint =  "Removing $Xpercent percent of the highest values ($nbValueToRemove values) gives you $resu as the longest intron. It's a good choice for MAKER ;-) \n";

print $ostreamReport $stringPrint;
if($opt_output){print $stringPrint;}

# Part 4
#########
# PLOT  #
#Choose a title:
my $title="Intron distribution";

## check using R
my $R = Statistics::R->new() or die "Problem with R : $!\n";

#calculate the breaks
my $breaks_ok=int($biggest_value/$opt_breaks);

#R command
$R->run(qq`
      
      listValues1=as.matrix(read.table("$pathIntron", sep="\t", he=F))
      pdf("$outputPDF")
      hist1<-hist(listValues1,breaks=$breaks_ok)
      plot(hist1\$mids,hist1\$counts)
      mylims <- par("usr")`
);

# Close the bridge
$R->stopR();



# remove temporary files
unlink $pathIntron;

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



__END__


=head1 NAME
 
gff3_sp_manage_introns.pl - This script give some information about introns (longest, shortest size mean ...) using the statistic method, 
then plot all the intron size values to get an overview of the introns size distribution.
It gives you as well the value of the longest intron after removing X percent(s) of the longest (removing potential biais / false positive).

=head1 SYNOPSIS

    ./gff3_sp_manage_introns.pl --gff=infile --out=outFile 
    ./gff3_sp_manage_introns.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>, B<-f>, B<--ref> or B<-reffile>

Input GFF3 file correponding to gene build. You can use several input files by doing: -f file1 -f file2 -f file3

=item  B<-w>, B<--window>, B<--break>, B<--breaks> or B<-b>

It the number of break used within the histogram plot. By default it's 1000. You can modify the value to get something more or less precise.

=item  B<-x>, B<--p>

Allows to modify the X values to calculate the percentage of the longest introns to remove. By default the value is 1 (We remove 1 percent of the longest).

=item  B<--out>, B<--output> or B<-o>

Output gff3 file where the gene incriminated will be write.

=item B<--help> or B<-h>

Display this helpful text.

=back

=cut
