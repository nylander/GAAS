#!/usr/bin/env perl

use warnings;
use strict;
use Pod::Usage;
use Bio::SeqIO ;
use Getopt::Long;
use File::Basename;
use GAAS::GAAS;

my $header = get_gaas_header();
my $infile = undef;
my $outfile = undef;
my $contigcount=0;
my $help = undef;

if ( !GetOptions( 'i=s' => \$infile,
                  'o|out|output=s' => \$outfile,
                  'h|help!'         => \$help ) )
{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($help) {
    pod2usage( { -verbose => 99,
                 -exitval => 0,
                 -message => "$header\n" } );
}

if ( ! defined( $infile) ) {
    pod2usage( {
           -message => "$header\nMust specify at least 1 parameters:\nInput fasta file (-i)\n",
           -verbose => 0,
           -exitval => 1 } );
}


#DEAL with output create default output file based on the input file name
if (! $outfile){
	 my ($file_in,$path,$ext) = fileparse($infile,qr/\.[^.]*/);
	 $outfile = $file_in.".agp";
}
open (AGP_FILE, ">$outfile");

# ===================

my $inseq = Bio::SeqIO->new('-file' => "<$infile",
               '-format' => 'Fasta' );

my $outseq = Bio::SeqIO->new(
            -file     => ">contigs.fasta",
            -format => 'fasta',
            );

#Read scaffolded FASTA-file
while (my $seq_obj = $inseq->next_seq ) {
  my $scaffold = $seq_obj->id;
  my $sequence = $seq_obj->seq;
  my $start=1;
  my $oldsum;
  my $newsum;
  my $count=0;
  my $rounded;

  next if ($scaffold =~ /^contig/i);
  foreach my $substring_sequence (split /(N{20,})/i, $sequence){
    my $type;
    my $substring_length = length($substring_sequence);
    $count++;
    $oldsum=$start;
    $newsum=$oldsum+$substring_length-1;

    if ($substring_sequence !~ m/^N+$/i){
      $type="W";
      $contigcount++;
 	  $rounded=sprintf("%05s", $contigcount);
 	  my $contig_obj = Bio::Seq->new(-seq => "$substring_sequence",
                                           -display_id => "contig$rounded",
                                           -alphabet => "dna" );
      $outseq->write_seq($contig_obj);
    }
    elsif ($substring_sequence =~ m/^N+$/i){
      $type="N";
    }
    $start += $substring_length;
    if ($type eq "W"){
     print AGP_FILE "$scaffold\t$oldsum\t$newsum\t$count\t$type\tcontig$rounded\t1\t$substring_length\t+\n";
    }
    if ($type eq "N"){
     print AGP_FILE "$scaffold\t$oldsum\t$newsum\t$count\t$type\t$substring_length\tscaffold\tyes\tpaired-ends\n";
    }
  }
}

close AGP_FILE;

__END__


=head1 NAME

gaas_scaffold2AGP.pl - This script

=head1 DESCRIPTION

Creates a AGP-file needed by e.g. EMBL for a scaffolded assembly

=head1 SYNOPSIS

    gaas_scaffold2AGP.pl -i scaffoldfile.fasta -o scaffoldfile.agp
    gaas_scaffold2AGP.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>, B<-f>, B<--ref> or B<-reffile>

Input fasta file.

=item  B<--out>, B<--output> or B<-o>

Output agp file.

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

AUTHOR - Jacques Dainat / Henrik Lantz
