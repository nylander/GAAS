#!/usr/bin/env perl

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use GAAS::GAAS;

my $header = get_gaas_header();
my $start_run = time();

my $inputFile;
my $output_suffix1=1;
my $output_suffix2=2;
my $gzip_output;
my $gzip_input;
my $pigz_compression_threads=1;
my $opt_help = 0;

Getopt::Long::Configure ('bundling');
if ( !GetOptions (
      'i|file|input|gff=s' => \$inputFile,
      'os1|output_suffix1=s' => \$output_suffix1,
      'os2|output_suffix2=s' => \$output_suffix2,
      'thread=i' => \$pigz_compression_threads,
      'c|gzip!' => \$gzip_output,
      'h|help!'         => \$opt_help )  )
{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

if ($opt_help) {
    pod2usage( { -verbose => 99,
                 -exitval => 0 } );
}

if (! $inputFile ){
   pod2usage( { -message => 'At least 1 input file is mandatory',
                 -verbose => 1,
                 -exitval => 1 } );
}

#Deal with extensions
my $pieces_list = split_keep_delimiter($inputFile);
my @pieces = @$pieces_list;
my $suffix  = pop(@pieces);
my $fq_ext;
my $filename;

if ($suffix eq ".gzip" or $suffix eq ".gz") {
  $gzip_input=1;
  $fq_ext = pop(@pieces);
  $filename = concat_list_from_left(\@pieces);
}
else{
  $fq_ext = $suffix;
  $filename = concat_list_from_left(\@pieces);
}

if ($gzip_input) {#unzip input case
	if ($gzip_output) {
    print "command1\n";
    my $command = 'gzip -dc '.$inputFile.' | paste - - - - - - - -  | tee >(cut -f 1-4 | tr "\t" "\n" | pigz --best --processes '.$pigz_compression_threads.' > '.$filename."_".$output_suffix1.$fq_ext.'.gz) | cut -f 5-8 | tr "\t" "\n" | pigz --best --processes '.$pigz_compression_threads.' > '.$filename."_".$output_suffix2.$fq_ext.".gz";
    print "Command launched:\n".$command."\n";
    system ("/bin/bash -c '$command'");
	}
	else{
    print "command2\n";
    my $command = 'gzip -dc '.$inputFile.' | paste - - - - - - - -  | tee >(cut -f 1-4 | tr "\t" "\n" > '.$filename."_".$output_suffix1.$fq_ext.') | cut -f 5-8 | tr "\t" "\n" > '.$filename."_".$output_suffix2.$fq_ext;
    print "Command launched:\n".$command."\n";
    system ("/bin/bash -c '$command'");
	}
}
else{
	if ($gzip_output) {
          print "command3\n";
          my $command = 'cat '.$inputFile.' | paste - - - - - - - -  | tee >(cut -f 1-4 | tr "\t" "\n" | pigz --best --processes '.$pigz_compression_threads.' > '.$filename."_".$output_suffix1.$fq_ext.'.gz) | cut -f 5-8 | tr "\t" "\n" | pigz --best --processes '.$pigz_compression_threads.' > '.$filename."_".$output_suffix2.$fq_ext.".gz";
          print "Command launched:\n".$command."\n";
          system ("/bin/bash -c '$command'");
        }
        else{
          print "command4\n";
          my $command = 'cat '.$inputFile.' | paste - - - - - - - -  | tee >(cut -f 1-4 | tr "\t" "\n" > '.$filename."_".$output_suffix1.$fq_ext.') | cut -f 5-8 | tr "\t" "\n" > '.$filename."_".$output_suffix2.$fq_ext;
          print "Command launched:\n".$command."\n";
          system ("/bin/bash -c '$command'");
        }
}

my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job done in $run_time seconds\n";


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

=head1 NAME

deinterleave_fastq.pl

=head1 SYNOPSIS

Deinterleaves a (compressed or not compressed) FASTQ file of paired reads into two FASTQ files.
Optionally GZip compresses the output FASTQ files using pigz.

Can deinterleave 100 million paired reads (200 million total
reads; a 43Gbyte file), in memory (/dev/shm), in 4m15s (255s)

Script inspired by a pure bash code from the nathanhaigh repository: https://gist.github.com/3521724
Also see the interleaving script: https://gist.github.com/4544979

The nathanhaigh script was itself inspired by Torsten Seemann's blog post:
http://thegenomefactory.blogspot.com.au/2012/05/cool-use-of-unix-paste-with-ngs.html

    deinterleave_fastq.pl -i input.fastq
    deinterleave_fastq.pl -i input.fastq --os R1 --os R2
    deinterleave_fastq.pl --help

The first command will create input.1.fastq and input.2.fastq files.
The second command will create input.R1.fastq and input.R2.fastq files.

=head1 OPTIONS

=over 8

=item B<-i>, B<--file> or B<--input>

STRING: Input fastq file that will be read.

=item B<--os1>, B<--output_suffix1>

STRING: Suffix to add to the output file 1. By default 1

=item B<--os2>, B<--output_suffix2>

STRING: Suffix to add to the output file 2. By default 2.

=item B<-gz> or B<--gzip>

Bolean: The output will be compressed using pigz.

=item B<--thread>

Integer: The number of thread used when running pigz.

=item B<--help> or B<-h>

Display this helpful text.

=back

=head1 FEEDBACK

=head2 Did you find a bug?

Do not hesitate to report bugs to help us keep track of the bugs and their
resolution. Please use the GitHub issue tracking system available at this
address:

            https://github.com/NBISweden/GAAS/issues

 Ensure that the bug was not already reported by searching under Issues.
 If you're unable to find an (open) issue addressing the problem, open a new one.
 Try as much as possible to include in the issue when relevant:
 - a clear description,
 - as much relevant information as possible,
 - the command used,
 - a data sample,
 - an explanation of the expected behaviour that is not occurring.

=head2 Do you want to contribute?

You are very welcome, visit this address for the Contributing guidelines:
https://github.com/NBISweden/GAAS/blob/master/CONTRIBUTING.md

=cut

AUTHOR - Jacques Dainat
