#!/usr/bin/env perl

###
# Implement case insensitive
###
use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Bio::SeqIO ;
use Bio::DB::Fasta;
use IO::File;
use GAAS::GAAS;

my $header = get_gaas_header();
my $start_run = time();
my $opt_fastafile;
my $opt_help = 0;
my $opt_OGfile;
my $opt_dir;
my $name_OG;
my %OG_seqID=();
my %OG_count=();
my @tab_seqID;
my $path;

# OPTION MANAGMENT
my @copyARGV=@ARGV;
if ( !GetOptions( 'f|fa|fasta=s' => \$opt_fastafile,
                  'og|OG_file=s' => \$opt_OGfile,
                  'd|dir=s' => \$opt_dir,
                  'h|help!'         => \$opt_help ) )
{
    pod2usage( { -message => "$header\nFailed to parse command line",
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($opt_help) {
    pod2usage( { -verbose => 99,
                 -exitval => 0,
                 -message => "$header\n" } );
}

if ( (! (defined($opt_OGfile)) ) or (! (defined($opt_fastafile)) ) ){
    pod2usage( {
           -message => "\nAt least 2 parametes are mandatory:\n Input reference fasta file (-f); orthoMCL group in text format obtained by the in-house orthoMCL pipeline (-og)\n\n".
           "Output will be created for you, one fasta file per OG group.\n",
           -verbose => 0,
           -exitval => 2 } );
}


##### MAIN ####

######################
####create folder
if (defined($opt_dir)){
  $path = "$opt_dir/";
}else {
  $path = "OG_fasta/";
}


my $createdir="mkdir -p $path";
if (-d $path)
  { print "$path already exists\n"; exit; }
else
  { print "mkdir $path \n" unless system($createdir) }


#### read fasta file and save info in memory
######################
my $db = Bio::DB::Fasta->new($opt_fastafile);
print ("Genome fasta parsed\n");

###########################
### open OG file

my $file = IO::File->new();
if ( defined $opt_OGfile ) {
  $file->open($opt_OGfile, 'r') or croak (sprintf("Can not open '%s' for reading: %s", $opt_OGfile, $!));
}

my $output_stat = "$path"."stat_OG.txt";

my $ostream_stat = IO::File->new();
$ostream_stat->open($output_stat, 'w' ) or
croak(
  sprintf( "Can not open '%s' for reading: %s", $output_stat, $! ) );

##############################
#create the hash with OG ID and the sequence IDs in the same orthology group
while(my $line =<$file>) {
  chomp $line;
  #print $line."\n";
  @tab_seqID = split (" ",$line);

  $tab_seqID[0]=~/^(OG_\d+)\:/;
  #print scalar @tab_seqID."\n";
  $name_OG=$1;

#chomp the first element of the table that is the OG ID
  @tab_seqID = @tab_seqID[ 1 .. $#tab_seqID ];

  $OG_seqID{$name_OG}=[@tab_seqID];

###to count the number of sequence per OG and species
  foreach (@tab_seqID) {
    $_=~/(\d+)\|/;
    $OG_count{$name_OG}{$1}++;
  }

}


##########################
#Now extract the sequences

foreach my $OG_ID (sort keys %OG_seqID){

  my $ostream_OG = "$path"."$OG_ID".".fasta";


  my $seqout_OG = Bio::SeqIO->new( -file   => ">$ostream_OG",
                                    -format => 'Fasta',
                                    );
  ####print number of sequences per OG

  print $ostream_stat "\n".$OG_ID."\t number of sequences:".scalar(@{$OG_seqID{$OG_ID}})."\n";
 #print $OG_ID."\t number of sequences: ".scalar(@{$OG_seqID{$OG_ID}})."\n";

    foreach my $count ( keys %{$OG_count{$OG_ID}}){
      ####print number of sequences per species and per OG
      print $ostream_stat $OG_ID."\t species: ".$count."\tnumber of sequences: ".$OG_count{$OG_ID}{$count}."\n";
  #  print $OG_ID."\t species: ".$count."\tnumber of sequences: ".$OG_count{$OG_ID}{$count}."\n";
    }

  foreach my $seq_name (@{$OG_seqID{$OG_ID}}) {

   if($db->seq($seq_name)){
    #create sequence object
    my $seq_obj = $db->get_Seq_by_id($seq_name);
    $seqout_OG->write_seq($seq_obj);

   }else{
     print $seq_name." not found into the $opt_fastafile fasta file !\n";
   }

 }


}

#END
print "usage: $0 @copyARGV\n";
my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job done in $run_time seconds\n";


__END__

=head1 NAME

gaas_fasta_extract_sequence_from_OG.pl

=head1 DESCRIPTION

This script extracts sequence in fasta format from a fasta file. You can extract one fasta sequence providing the name of a file created by the in-house orthoMCL pipeline.
The OG file contains all the orthoMCL groups and the ID of the sequences in each group.

=head1 SYNOPSIS

    gaas_fasta_extract_sequence_from_OG.pl -f infile.fasta -og OGfile.txt
    gaas_fasta_extract_sequence_from_OG.pl --help

=head1 OPTIONS

=over 8

=item B<-f> or B<--fasta>

Input fasta file.

=item B<-og>, B<--og>

The OG file contains all the orthoMCL groups and the ID of the sequences in each group.

eg :

OG_1000: 5833|MAL13P1.2:pep 5833|PF10_0398:pep

OG_1001: 5833|MAL13P1.1:pep 5833|PFE0005w:pep 5833|MAL8P1.220:pep 5833|PFF1595c:pep

=item B<-d>, B<--dir>

optional you can choose a name for the output folder, by default it will be called OG_fasta

=item B<-h> or B<--help>

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
