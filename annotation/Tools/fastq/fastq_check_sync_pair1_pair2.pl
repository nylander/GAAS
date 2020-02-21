#!/usr/bin/env perl

###################################################
# Jacques Dainat 01/2018                          #
# jacques.dainat@nbis.se                          #
###################################################

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Data::Dumper;

my $start_run = time();

my @inputFile;
my $check_complete;
my $gzip_input;
my $nb;
my $opt_help = 0;

Getopt::Long::Configure ('bundling');
if ( !GetOptions (
      'i|file|input=s'  => \@inputFile,
      'c|complete!'     => \$check_complete,
      'nb=s'            => \$nb,
      'h|help!'         => \$opt_help )  )
{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

if ($opt_help) {
    pod2usage( { -verbose => 2,
                 -exitval => 0 } );
}

if (! ($#inputFile == 1) ){
   pod2usage( { -message => 'at least 2 input files are mandatory',
                 -verbose => 1,
                 -exitval => 1 } );
}

#Deal with extensions
#my $pieces_list = split_keep_delimiter($inputFile);
#my @pieces = @$pieces_list;
#my $suffix  = pop(@pieces);

#if ($suffix eq ".gzip" or $suffix eq ".gz") {
#  $gzip_input=1;
#}

my $in1 = IO::File->new();
my $in2 = IO::File->new();
open($in1, '<', $inputFile[0]) or die "Could not open file '$inputFile[0]' $!";
open($in2, '<', $inputFile[1]) or die "Could not open file '$inputFile[1]' $!";


my $read_cpt = 0;
my $read_fail = 0;
my $header_type=undef;
my $count=0;

while (!eof($in1) and !eof($in2)) {
  $count++;
  my $id1 = <$in1>;
  my $id2 = <$in2>;
  # skip all line that are not header 
  next unless ($count % 4 == 1 ); 

  #extract header
  chomp $id1;
  chomp $id2;

  #check header type
  if ($read_cpt == 0){
    if ($id1 =~ /^@\S+\/[12]$/) { # @ at the start of the line followed by non-whitespace, a /, a 1 or 2, the end of the line
      $header_type = 1;
      print STDOUT "Read Id looks like Casava 1.7 style\n"; # TESTING
    } 
    elsif ($id1 =~ /^@\S+\W[12]\S+$/) { # @ at the start of the line followed by non-whitspace, a space, a 1 or 2, non-whitespace
      $header_type = 2;
      print STDOUT "Read Id looks like Casava 1.8 style\n"; 
    } 
    else {
      print STDOUT "Unknwon id style (Not Casava 1.7 or 1.8): $id1\n";
      exit 1;
    }
  }


  if ($header_type == 1) { # 1 or 2 at end of id
      chop $id1;   # last char of the id (should be the "1" or "2")
      chop $id2;
  } else {
      ($id1) = split ' ', $id1;
      ($id2) = split ' ', $id2;
  }

  if ($id1 ne $id2 ){
    if($check_complete) {
      $read_fail++;
    }
    else{
      my $end_run = time();
      my $run_time = $end_run - $start_run;
      print "ERROR dectected. Read1 and Read2 not synchronized at line $count.\nRead1 = $id1 -- Read2 = $id2\n";
      print "Runtime $run_time\n";
      exit;
    }
  }
  $read_cpt++;

  if($nb and $nb == $read_cpt){
    last;
  }
}

my $percent = ($read_fail / $read_cpt) * 100 ;
$percent = $percent * 100; # make it percent
$percent = sprintf("%.2f", $percent);

my $end_run = time();
my $run_time = $end_run - $start_run;

print "Check successfully passed in $run_time seconds. $read_cpt reads read.\n";

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


sub concat_list_from_left{

  my ($list)=@_;

  my $result="";
  foreach my $element (@{$list}){
    $result = $result.$element;  
  }

  return $result;
}

sub split_keep_delimiter{
  my ($string)=@_;

  my @result;
  my @pieces = split(/\./, $string);
  my $cpt=0;
  foreach my $element (@pieces){
    if($cpt != 0){
      push @result, ".".$element;
    }
    else{
      push @result, $element;
    }
    $cpt++;
  } 
  return \@result;
}

__END__

      'c|complete!' => \$check_complete,
      'h|help!'         => \$opt_help )  

=head1 NAME

fastq_check_sync_pair1_pair2.pl

=head1 SYNOPSIS

The aim of this script is to check that paired reads from 2 fastq files are still synchronized. 
Read1 and the read2, that come from a paired sequencing, are in the same position in the two fastq files. 
But the order of read in R1/R2 files can get out of sync if you e.g scan/trim the two files independently. So it is a good thing always to check.

    fastq_check_sync_pair1_pair2.pl -i input_R1.fastq -i input_R2.fastq
    fastq_check_sync_pair1_pair2.pl --help

=head1 OPTIONS

=over 8

=item B<-i>, B<--file> or B<--input>

STRING: Input fastq file that will be read.

=item B<--complete> or B<-c>

BOLEAN - In complete mode, the script doesn't stop at the first synchronization problem, but will read the whole file and report the number of de-synchronization found.

=item B<--nb>

Integer - Allow to check just a subsample of the reads. So, define here the number of read to check.

=item B<--help> or B<-h>

Display this helpful text.

=back

=cut
