#!/usr/bin/env perl


use Carp;
use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO;
use GAAS::GAAS;

my $header = get_gaas_header();
my $outfile = undef;
my $file1 = undef;
my $file2 = undef;
my $help= 0;

if ( !GetOptions(
    "help|h" => \$help,
    "fasta1|file1|f1=s" => \$file1,
    "fasta2|file2|f2=s" => \$file2,
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

if ( ! ((defined($file1)) and (defined($file2)))){
    pod2usage( {
           -message => "$header\nAt least 2 parameters are mandatory.\n",
           -verbose => 0,
           -exitval => 1 } );
}

######################
# Manage output file #
my $fastaout;
if ($outfile) {
  $outfile=~ s/.gff//g;
open(my $fh, '>', $outfile.".gff") or die "Could not open file '$outfile' $!";
  $fastaout=  Bio::SeqIO->new(-fh => $fh , '-format' => 'Fasta');
}
else{
  $fastaout = Bio::SeqIO->new(-fh => \*STDOUT, -format => 'Fasta');
}


                #####################
                #     MAIN          #
                #####################


########################
### Manage INPUT FILES #
my $fasta1  = Bio::SeqIO->new(-file => $file1 , -format => 'Fasta');
my $fasta2  = Bio::SeqIO->new(-file => $file2 , -format => 'Fasta');

##############
#### MAIN ####
my $nbToRemove;
my $nbRemoved;
my $nbSeqPrint;
my %id_fasta1;
while ( my $seq = $fasta1->next_seq() ) {
 $id_fasta1{$seq->id}++;
 $nbToRemove++;
}

while ( my $seq = $fasta2->next_seq() ) {
 if(! exists($id_fasta1{$seq->id})){
    $fastaout->write_seq($seq);
    $nbSeqPrint++;
 }
 else{$nbRemoved++;}
}

my $totalSeq=$nbRemoved+$nbSeqPrint;
print "On the $nbToRemove sequences in $file1, $nbRemoved sequences have been removed from $file2.\nSo, on the $totalSeq sequences of $file2, $nbSeqPrint have been printed.\n";

__END__

=head1 NAME

Compare two fasta file in order to remove occurence of fasta sequence from file 1 present in file 2.
The whole header must be identical to be consider as identic.

=head1 SYNOPSIS

    perl my_script.pl --fasta1 file1 --fasta2 file2 [--out outfile]
    perl my_script.pl --help

=head1 OPTIONS

=over 8

=item B<--fasta1>, B<--file1> or B<-f1>

Fasta file 1. The headers of sequences of this file will be used to compare against those to file 2.

=item B<--fasta2>, B<--file2> or B<-f2>

Fasta file 2. This is the "reference file" in which we will remove sequences already existing in file 1.

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
