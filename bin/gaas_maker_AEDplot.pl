#!/usr/bin/env perl

#############################################
# Jacques Dainat 11/2014 # This version use the GFF.pm from Andreas code
#############################################

#libraries
use Statistics::R;
use File::Basename;
use strict;
use warnings;
use Data::Dumper;
use Carp;
use Getopt::Long;
use IO::File;
use Pod::Usage;
use List::MoreUtils qw(uniq);
use Bio::Tools::GFF;
use GAAS::CheckModule qw(:Ok);
# END libraries

# PARAMETERS - OPTION
my @opt_files;
my $opt_output;
my $opt_breaks;
my $opt_help = 0;
# END PARAMETERS - OPTION


# OPTION MANAGMENT
if ( !GetOptions( 'f=s' => \@opt_files,
                  'w|window=i'      => \$opt_breaks,
                  'o|output=s'      => \$opt_output,
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

if ( ! ($#opt_files >= 0)){
    pod2usage( {
           -message => "\nAt least 1 parameter is mandatory:\nInput gff file\n\n".
           "You may add as many file you want like: -f file1 -f file2 -f file3\n".
           "Many optional parameters are available. Look at the help documentation to know more.\n",
           -verbose => 0,
           -exitval => 2 } );
}

#############################
####### Manage options #######
#############################

#Choose breaks value:
if(! $opt_breaks){
  $opt_breaks="0.05";
}

#############################
####### Manage output #######
#############################
my $outputPDF;
if (defined($opt_output) ) {
  if (-f $opt_output){
      print "Cannot create a directory with the name $opt_output because a file with this name already exists.\n";exit();
  }
 $outputPDF=$opt_output."pdf";
}
else{
  $outputPDF="outputPlot.pdf";
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


#######################
#        MAIN         #
#######################

########## constant #############
my $R_command;
my $nbFile=0;
my @listTmpFile;
#Choose a title:
my $title="AED distribution";

#PART 1
###################################
# Read input gff3 files one by one and save value in hash of list
my %hashOfList;
foreach my $file (@opt_files){

  my $gffio = Bio::Tools::GFF->new(-file => $file, -gff_version => 3);

  # parse file name te remove extension
  my ($file1,$dir1,$ext1) = fileparse($file, qr/\.[^.]*/);

  #Parse GFF to get AED information for each mRNA
  my $listRef=parseGFF($gffio);

  $hashOfList{$file1}=$listRef;
  $gffio->close();

  print("Parsing $file Finished\n\n");
}

#PART 2
###############################
#print values in file and pre-plot to get highest Y axis.
my $highestYaxis=0;
foreach my $fileName (keys %hashOfList){

  #########################
  # Write value in tmp files
  my $pathAED="tmp_AED_".$fileName.".txt";
  push (@listTmpFile, $pathAED);
  # Manage Output
  my $ostreamAED   = IO::File->new();
  $ostreamAED->open( $pathAED, 'w' ) or
        croak(
          sprintf( "Can not open '%s' for writing %s", $pathAED, $! )
        );
  foreach  my $AEDvalue (@{$hashOfList{$fileName}}){
    print $ostreamAED "$AEDvalue\n";
  }
  $ostreamAED->close();

  ## check using R
  my $R = Statistics::R->new() or die "Problem with R : $!\n";

  #R command
  $R->run(qq`

        listValues1=as.matrix(read.table("$pathAED", sep="\t", he=F))

        # create a break point list formated correctly in purpose
        a<-seq(0,0.9999,$opt_breaks)
        a[length(a)+1]<-0.99999
        a[length(a)+1]<-1
        breakingPointList<-c(0,a)

        hist1<-hist(listValues1, breaks=breakingPointList, plot=F)
        plot(hist1\$mids,hist1\$counts)
        #par(new=TRUE)
        mylims <- par("usr")`
  );

  #retrieve R values in Perl
  my $maxY = $R->get('mylims');

  if($maxY->[$#$maxY] > $highestYaxis){
    $highestYaxis=$maxY->[$#$maxY];
  }

  # Close the bridge
  $R->stopR();
}

#PART 3
############################################
##  Read values from files and plot with correct Y axis
foreach my $fileName (keys %hashOfList){
  $nbFile++;
  my $pathAED="tmp_AED_".$fileName.".txt"; # Need to be similar as in #PART 2


  ##################################
  # create main part of R command  #
  ##################################
  if($nbFile > 1){ # only one plot to do
    $R_command.='
                par(new=TRUE)';
                # write it for each file
  $R_command.=write_R_command($pathAED, $fileName, $highestYaxis, $nbFile, $outputPDF);
  }
  else{
    # write it for each file
    $R_command.=write_first_R_command($pathAED, $fileName, $highestYaxis, $nbFile, $outputPDF);
  }
}


####################################
# Surround main part of R command  #
####################################

  #add header
  my $final_R_command='#create output
                  pdf("'.$outputPDF.'")
                  # create an empty vector
                  listlegend <- c();';

  #add heart
  $final_R_command.=$R_command;

  #add footer
  $final_R_command .='# Add Title
                title(main="'.$title.'")

                #Add Legend
                legend("topright", col=(1:'.$nbFile.'), lty=1, c(listlegend))';

# plot
plotR($final_R_command);

# remove temporary files
unlink @listTmpFile;

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

sub write_first_R_command{
  my ($pathIn1,$name1,$yAxisValue,$colorNb,$outputPDF)=@_;

  my $command='

      listValues1=as.matrix(read.table("'.$pathIn1.'", sep="\t", he=F))
      legendInfo=paste("'.$name1.'","(",length(listValues1),"mRNAs )")
      listlegend<-c(listlegend,legendInfo)

      # create a break point list formated correctly in purpose
      a<-seq(0,0.9999,'.$opt_breaks.')
      a[length(a)+1]<-0.99999
      a[length(a)+1]<-1
      breakingPointList<-c(0,a)

      hist1<-hist(listValues1, breaks=breakingPointList, plot=F)
      plot(hist1$mids,hist1$counts, type="l", ylim=c(0,'.$yAxisValue.'), col='.$colorNb.', main="", xlab="AED score", ylab="Number of mRNA")
      ';

      return $command
}

#xlab="AED score"
sub write_R_command{
  my ($pathIn1,$name1,$yAxisValue,$colorNb,$outputPDF)=@_;

  my $command='

      listValues1=as.matrix(read.table("'.$pathIn1.'", sep="\t", he=F))
      legendInfo=paste("'.$name1.'","(",length(listValues1),"mRNAs )")
      listlegend<-c(listlegend,legendInfo)

      # create a break point list formated correctly in purpose
      a<-seq(0,0.9999,'.$opt_breaks.')
      a[length(a)+1]<-0.99999
      a[length(a)+1]<-1
      breakingPointList<-c(0,a)

      hist1<-hist(listValues1, breaks=breakingPointList, plot=F)
      plot(hist1$mids,hist1$counts, type="l", ylim=c(0,'.$yAxisValue.'), col='.$colorNb.', main="", yaxt="n", xaxt="n", xlab="", ylab="")
      ';

      return $command
}

sub plotR {
  my ($command)=@_;

  my $R = Statistics::R->new() or die "Problem with R : $!\n";

#R command
$R->send(
      qq`

      $command

      dev.off()`
          );

# Close the bridge
$R->stopR();
}


# Delete temporary file
#unlink "$pathIn1";
#unlink "$pathIn2";

# method to parse GFF3 files
# take in account features gens,mRNA,tRNA,exon,CDS,three_prime_UTR and five_prime_UTR
sub parseGFF {
  my @list;
  my($file_in) = @_;
  print( "Reading features from $file_in...\n");
  # read file and decompose it
  while (my $feature = $file_in->next_feature() ) {

        my $type = $feature->primary_tag();

        if (lc($type) eq 'mrna'){

          if(! $feature->has_tag('_AED')){
              print "AED of this feature not found for".$feature->_tag_value('ID')."\n";
          }
          else{
             my $AED=$feature->_tag_value('_AED');
             push(@list,$AED);
          }
        }
  }
  return \@list;
}

__END__

=head1 NAME

AEDplot.pl -
The script take one or several gff file(s) as input from Maker and create a Plot of their AED score (Attributes used: "_AED"). -
=head1 SYNOPSIS

    ./maker_AEDplot.pl -f infile1.gff[ --output outfile ]
    ./maker_AEDplot.pl --help

=head1 OPTIONS

=over 8

=item B<-f>

Input GFF3 file(s) created with maker (with mRNA containing _AED attribute). When you want to use several distinct files do: -f file1 -f file2 -f file3

=item B<-w>, B<--window>
The AED score value is between 0 and 1. You can define how precise the plot will by defining the window size that will be taken in account to peform the calcul. The value is 0.05 By default.

=item B<--output>, B<-o>

Output name of the pdf file created. If none provided, the default output is ouputPlot.pdf

=item B<-h> or B<--help>

Display this helpful text.

=back

=cut
