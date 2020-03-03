#!/usr/bin/env perl

#LIBRARIES
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use IO::File;
use GAAS::GAAS;

my $header = get_gaas_header();
#PARAMETERS
my $input_listcontigWrong;
my $input_datastorelog;
my $output_file;
my $delete_contig;
my %contigs;
my %datastores;

#get arg and PARAMETERS
{
  my ($input_listcontigWrong, $input_datastorelog, $output_file, $delete_contig);

  # Define script options
  GetOptions(
  'help|h'                => sub { pod2usage( -verbose => 99 , -exitval => 0, -message => "$header\n" )},
  'output|o=s'       => \$output_file,
  'contig-name|f=s'       => \$input_listcontigWrong,
  'datastore-name|d=s'    => \$input_datastorelog,
  'delete-contig|c=s'    => \$delete_contig,

  ) or pod2usage(2);

  pod2usage( "--contig-name must be specified" )
  unless defined $input_listcontigWrong;

  pod2usage( "--datastore-name must be specified" )
  unless defined $input_datastorelog;

  pod2usage( "--delete-contig options must be specified" )
  unless defined $delete_contig;

  main($input_listcontigWrong, $input_datastorelog, $output_file, $delete_contig);
}

sub main {
  my ($input_listcontigWrong, $input_datastorelog, $output_file, $delete_contig) = @_;

  my $ostream = IO::File->new();
  if(defined($output_file)){
    $ostream->open($output_file, 'w' ) or
    croak(
      sprintf( "Can not open '%s' for reading: %s", $output_file, $! ) );
    }
  else{
    $ostream->fdopen( fileno(STDOUT), 'w' ) or
      croak( sprintf( "Can not open STDOUT for writing: %s", $! ) );
  }
get_list($input_listcontigWrong, $input_datastorelog, $ostream, $delete_contig);

}

sub get_list {

  my ($input_listcontigWrong, $input_datastorelog, $ostream, $delete_contig) = @_;

  my $listcontigWrong = IO::File->new();
  if ( defined $input_listcontigWrong ) {
    $listcontigWrong->open($input_listcontigWrong, 'r') or croak (sprintf("Can not open '%s' for reading: %s", $input_listcontigWrong, $!));
  }

  my $datastorelog = IO::File->new();
  if ( defined $datastorelog ) {
    $datastorelog->open($input_datastorelog, 'r') or croak (sprintf("Can not open '%s' for reading: %s", $input_datastorelog, $!));
  }

  #get list of path where contigs are and delete them
  while(my $line =<$listcontigWrong>) {
    chomp $line;
    if ($line=~/(genome.maker.output\/)(genome_datastore\/\S+\/)/){
      $contigs{$2}=$line;
      if ($delete_contig eq 'all') {
        system("rm -r $1$2");
        }
      }
    }
  #then recreate a file without lines corresponding to the contigs in output of maker (genome_master_datastore_index.log)
  while(my $lineB =<$datastorelog>) {
    chomp $lineB;
    if ($lineB=~/(genome_datastore\/\S+\/)/){
      if ($delete_contig eq 'all' || $delete_contig eq 'log'){
        if(! exists $contigs{$1}) {
          print $ostream $lineB."\n";
        }
      }
    }
  }
}

1;

__END__;

=head1 NAME

gaas_maker_get_rid_of_contig.pl

=head1 DESCRIPTION

Get rid of contigs not processed properly by maker in the log file and in the output folders of maker.
Create a new log file that will need to be renamed as genome_master_datastore_index.log to replace the old one.
Then maker can be rerun with this new log file.

=head1 SYNOPSIS

Get rid of contigs not processed properly by maker in the log file and in the output folders of maker.
Create a new log file that will need to be renamed as genome_master_datastore_index.log to replace the old one.
Then maker can be rerun with this new log file.

gaas_maker_get_rid_of_contig.pl --help

gaas_maker_get_rid_of_contig.pl --contig-name|-f

gaas_maker_get_rid_of_contig.pl --output-name|-o

gaas_maker_get_rid_of_contig.pl --datastore-name|-d

gaas_maker_get_rid_of_contig.pl --delete-contig|-c log/all

=head1 OPTIONS

=over 8

=item B<--datastore-name> or B<-d>

Input datastore log file

=item B<--contig-name> or B<-f>

Input file containing the list of wrong contig

=item B<--delete-contig> or B<-c>

<log> option will only delete contigs in the log file
<all> option will delete contigs in the log file and contigs' folders

=item B<--output> or B<-o>

File output name

=item B<--help>

Display the manual page.

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

AUTHOR - Lucile Soler / Jacques Dainat
