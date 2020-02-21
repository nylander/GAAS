#!/usr/bin/env perl

#LIBRARIES
use strict;
use warnings;
use Data::Dumper;
use YAML::XS 'LoadFile';
use Getopt::Long;
use Pod::Usage;
use IO::File;

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
  'help|h'                => sub { pod2usage( -verbose => 2 ), -exitval => 0 },
  'man'                   => sub { pod2usage( -verbose => 2 )},
  'output-name|o=s'       => \$output_file,
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

Get rid of contigs not processed properly by maker in the log file and in the output folders of maker.
Create a new log file that will need to be renamed as genome_master_datastore_index.log to replace the old one.
Then maker can be rerun with this new log file.

=head1 AUTHOR

Lucile SOLER NBIS 18/07/2016

=head1 SYNOPSIS

Get rid of contigs not processed properly by maker in the log file and in the output folders of maker.
Create a new log file that will need to be renamed as genome_master_datastore_index.log to replace the old one.
Then maker can be rerun with this new log file.

getRidOfContig.pl --help

getRidOfContig.pl --contig-name|-f

getRidOfContig.pl --output-name|-o

getRidOfContig.pl --datastore-name|-d

getRidOfContig.pl --delete-contig|-c log/all

    log option will only delete contigs in the log file
    all option will delete contigs in the log file and contigs' folders

=head1 OPTIONS

=over

=item B<--help>

Display a brief usage message.

=item B<--man>

Display the manual page.

=back

=head1 DESCRIPTION

Get rid of contigs not processed properly by maker in the log file and in the output folders of maker.
Create a new log file that will need to be renamed as genome_master_datastore_index.log to replace the old one.
Then maker can be rerun with this new log file.

=cut
