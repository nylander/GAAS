#!/usr/bin/env perl

use Carp;
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO;
use IO::File;
use GAAS::GAAS;

my $header = get_gaas_header();
my $outfile = undef;
my $file1 = undef;
my $agp = undef;
my $verbose= undef;
my $help= 0;

if ( !GetOptions(
    "help|h" => \$help,
    "fasta|file|f=s" => \$file1,
    "a|agp=s" => \$agp,
    "output|outfile|out|o=s" => \$outfile))

{
    pod2usage( { -message => "$header"."Failed to parse command line.",
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($help) {
    pod2usage( { -message => "$header\n",
                 -verbose => 99,
                 -exitval => 0 } );
}

if ( ! ((defined($file1)) and (defined($agp)))){
    pod2usage( {
           -message => "$header\nAt least 2 parameters are mandatory.\n",
           -verbose => 0,
           -exitval => 1 } );
}

######################
# Manage output file #
my $fastaout;
if ($outfile) {
  $outfile=~ s/.fasta//g;
  $outfile=~ s/.fa//g;
open(my $fh, '>', $outfile.".fa") or die "Could not open file '$outfile' $!";
  $fastaout=  Bio::SeqIO->new(-fh => $fh , -format => 'Fasta');
}
else{
  $fastaout = Bio::SeqIO->new(-fh => \*STDOUT, -format => 'Fasta');
}


                #####################
                #     MAIN          #
                #####################

######################
### Parse AGP input #
my @list_contig_used;
print "Reading file $agp\n";
if (open(my $fh, '<:encoding(UTF-8)', $agp)) {
  while (my $row = <$fh>) {
    if($row =~ /^#/){next;}
    chomp $row;
    my ($object, $object_beg, $object_end, $part_number, $component_type, $component_id_or_gap_length,
        $component_beg_or_gap_type, $component_end_or_linkage, $orientation_or_linkage_evidence ) = split(/\t/, $row);
    if($component_type ne "N" and $component_type ne "U"){
      if ($object =~ /^Super-Scaffold/ ){
        #print  $component_id_or_gap_length."\n";
        push (@list_contig_used, $component_id_or_gap_length);
      }
    }
  }
}
else {
  warn "Could not open file $agp $!";
}

print "Parsing Finished\n";
### END Parse AGP input #
#########################

# primary contig header: >004069F|arrow|arrow_obj
# alternative contig header: >000019F-023-01|arrow|arrow_obj
# If piece used not all contig: >000838F|arrow|arrow_subseq_300884:420004 => To keep


#Go all over fasta1 and skip sequence to exclude
my $fasta1  = Bio::SeqIO->new(-file => $file1 , -format => 'Fasta');
my $nbRemoved=0;
while ( my $seq = $fasta1->next_seq() ) {
  my $header = $seq->id;
  if( $header =~ m/-/ ){ #its from alternative
    #print "alternative:  $header\n";
    if( $header =~ m/subseq/ ){
      #print "Subsequence we have to keep it\n";
      if ( $header =~ m/(.+?(?=-))/ ) {
        #print "header1= $header $1\n";
        if ( grep( /^\Q$1\E\|/, @list_contig_used ) ) {
          my @res = grep( /^\Q$1\E\|/, @list_contig_used );
          #print "match in AGP: @res\n";
          if(@res > 1 ){
            print "header used severeal times in the agp file!!! \n";
          }
          else{
            print "The primary version of this contig (@res) has been already taken into the hybrid_Assembly.fasta output. No need to include it \n" if($verbose);
          }
        }
        else{
          print "subseq of alternative contigs $header has to be included into the final assembly. Indeed the primary version of this one has not been included\n";
          $fastaout->write_seq($seq);
        }
      }
    }
  }
  else{ #its from primary
    if( $header =~ m/subseq/ ){
      print "Subsequence we have to keep it\n" if ($verbose);
      $fastaout->write_seq($seq);
    }
    else{
      if ( $header =~ m/(.+?(?=_obj))/ ) {
        #print "header= $header $1\n";
        if ( grep( /^\Q$1\E/, @list_contig_used ) ) {
          my @res = grep( /^\Q$1\E/, @list_contig_used );
          print "match in AGP: @res\n";
          if(@res > 1 ){
            print "header used severeal times in the agp file!!! Need to implement how to deal with this case. A loop will be enough...\n";
          }
          else{
            print "header used in AGP: @res \n";
          }
        }
        else{
           $fastaout->write_seq($seq);
        }
      }
      else{
        print "not match for <_obj> at the end of the string\n";
      }
    }
  }
}



__END__

=head1 NAME

gaas_fasta_bionano_filter_not_scaffolded_part.pl

=head1 DESCRIPTION

This script aims to filter the NOT_SCAFFOLDED.fasta file from bionano output in order to remove redundant part from secondary assembly. Indeed the NOT_SCAFFOLDED.fasta file is a mixup of the primary and the secondary assembly.

Is not included in the output:
  - piece of the secondary assembly (they are cut in pieces when a piece is used into a scaffold), when the corresponding sequence from the primary assembly is already used.
  - contig of the secondary assembly
  - contig of the primary assembly if the counterpart of the secondary assembly is already used into a scaffold.

=head1 SYNOPSIS

    gaas_fasta_bionano_filter_not_scaffolded_part.pl my_script.pl --fasta1 file1 -a agp [--out outfile]
    gaas_fasta_bionano_filter_not_scaffolded_part.pl my_script.pl --help

=head1 OPTIONS

=over 8

=item B<--fasta1>, B<--file1> or B<-f1>

Fasta file 1. The headers of sequences of this file will be used to compare against those to file 2.

=item B<-a>, B<--agp> or B<-f2>

This is a file containing the headers of sequence to be removed. Only one ID per line. Header should be identical at 100% to be removed.

=item B<-o> , B<--output> , B<--out> or B<--outfile>

Output fasta file.  If no output file is specified, the output will be
written to STDOUT.

=item B<--help> or B<-h>

Getting help.
Display the full information.

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
