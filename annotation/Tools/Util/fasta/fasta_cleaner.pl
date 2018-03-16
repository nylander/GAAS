#!/usr/bin/env perl

# A filter for Uniprot and RefSeq fasta files that makes the fasta
# headers a bit more terse.  Reads from STDIN, writes to STDOUT.
#
## Note: Will pass any other fasta file unchanged.
## Note: For Uniprot, will also change any 'O' in the protein sequence
##       into 'K'.

use strict;
use warnings;
use Pod::Usage;
use POSIX qw(strftime);
use Getopt::Long;

my $start_run = time();
my $file_fasta;
my $outfile;
my $help = 0;
my $verbose = 0;
my $header = qq{
########################################################
# BILS 2018 - Sweden                                   #  
# jacques.dainat\@nbis.se                               #
# Please cite NBIS (www.nbis.se) when using this tool. #
########################################################
};

my @copyARGV=@ARGV;
Getopt::Long::Configure ('bundling');

if ( !GetOptions(
    "help|h" => \$help,
    "fasta|fa|f=s" => \$file_fasta,
    "v!" => \$verbose,
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
 
if ( !(defined($file_fasta)) ){
    pod2usage( {
           -message => "$header\nAt least 1 parameter is mandatory:\n Input fasta file (--fasta)\n\n",
           -verbose => 0,
           -exitval => 1 } );
}

open my $fh, $file_fasta or die "Could not open $file_fasta: $!";

#OUTPUT
my $fho;
if ($outfile) {
  open($fho, '>', $outfile) or die "Could not open file '$outfile' $!";
  }
else{
  $fho = *STDOUT;
}

# To follow progression
  my $startP=time;
  my $nbLine=`grep -c ">" $file_fasta`;
  $nbLine =~ s/ //g;
  chomp $nbLine;
  print "$nbLine sequence to process...\n";
  my $line_cpt=0;

#########
#MAIN

my $first_line = <$fh>;
chomp($first_line);

if ( $first_line !~ /^>/ ) {
    die( sprintf( "This does not look like fasta formatted input:\n%s\n",
                  $first_line ) );
}

my %parsers = (
    'null' => sub {
        my ($line) = @_;
        return $line;
    },

    'uniprot' => sub {
        my ($line) = @_;


        if ( $line =~ /^>(?:sp|tr)\|([^|]+).*PE=(\d+) SV=(\d+)/ ) {
            $line_cpt++;
            return sprintf( ">%s.%d %d", $1, $3, $2 );
        }
        else {
            $line =~ tr/O/K/;
        }

        return $line;
    },
    'refseq' => sub {
        my ($line) = @_;

        if ( $line =~ /^>gi/ ) {
            $line_cpt++;
            return sprintf( ">%s", [ split( /\|/, $line ) ]->[3] );
        }

        return $line;
    } );

my $parser = 'null';

if    ( $first_line =~ /^>(?:sp|tr)/ ) { $parser = 'uniprot'; }
elsif ( $first_line =~ /^>gi/ ) { $parser = 'refseq'; }

print $fho ( $parsers{$parser}($first_line), "\n" );

while ( my $line = <$fh> ) {
    chomp($line);
    print $fho $parsers{$parser}($line), "\n";

    #Display progression
    if ((30 - (time - $startP)) < 0) {
      my $done = ($line_cpt*100)/$nbLine;
      $done = sprintf ('%.0f', $done);
          print "\rProgress : $done %";
      $startP= time;
    }
}

#END
print "usage: $0 @copyARGV\n";
my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job done in $run_time seconds\n";