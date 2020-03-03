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
my $file2 = undef;
my $help= 0;

if ( !GetOptions(
    "help|h" => \$help,
    "fasta|file|f=s" => \$file1,
    "list|l=s" => \$file2,
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
                 -exitval => 0
                } );
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


########################
### Manage INPUT FILES #
my $nbToExclude;
my %list_to_exclude;
if (-f $file2){
	my $ID_list  = IO::File->new("<".$file2);

	##############
	#### MAIN ####

	# create hash of ID to remove
	while ( <$ID_list> ) {
	  chomp;
	  if(! $_ =~ /^\s*$/){
	    $list_to_exclude{$_}++;
	    $nbToExclude++;
	  }
	}
}
else{
	$list_to_exclude{$file2}++;
        $nbToExclude++;
}
print "You want to removed $nbToExclude sequence from $file1\n";

#Go all over fasta1 and skip sequence to exclude
my $fasta1  = Bio::SeqIO->new(-file => $file1 , -format => 'Fasta');
my $nbRemoved=0;
while ( my $seq = $fasta1->next_seq() ) {
if(! exists($list_to_exclude{$seq->id})){
    $fastaout->write_seq($seq);
 }
 else{$nbRemoved++;}
}

if($nbToExclude == $nbRemoved){
  print "Exclusion successful. All protein you wanted to exclude have not been kept in the output.\n";
}else{
  print "WARNING only $nbRemoved sequences on the $nbToExclude you wanted to exclude have been excluded fron the output\n";
}


__END__

=head1 NAME

gaas_fasta_removeSeqFromIDlist

=head1 DESCRIPTION

Compare a fasta file to a list of ID in order to remove the matching name from file 1.
The whole header must be identical to be consider as identic.

=head1 SYNOPSIS

    gaas_fasta_removeSeqFromIDlist.pl --fasta1 file1 --list file2 [--out outfile]
    gaas_fasta_removeSeqFromIDlist.pl --help

=head1 OPTIONS

=over 8

=item B<--fasta1>, B<--file1> or B<-f1>

Fasta file 1. The headers of sequences of this file will be used to compare against those to file 2.

=item B<--fasta2>, B<--file2> or B<-f2>

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
