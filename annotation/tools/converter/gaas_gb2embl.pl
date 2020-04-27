#!/usr/bin/env perl

## NBIS 2015
## jacques.dainat@nbis.se

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Bio::SeqIO;
use GAAS::GAAS;

my $header = get_gaas_header();
my $outfile = undef;
my $gb = undef;
my $help;

if( !GetOptions(
    "help" => \$help,
    "gb=s" => \$gb,
    "outfile|output|o|out|embl=s" => \$outfile))
{
    pod2usage( { -message => "Failed to parse command line\n$header\n",
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($help) {
    pod2usage( { -verbose => 99,
                 -exitval => 0,
                 -message => "$header\n" } );
}

if ( ! (defined($gb)) ){
    pod2usage( {
           -message => "$header\nMissing the --gb argument\n",
           -verbose => 0,
           -exitval => 1 } );
}

## Manage output file
my $embl_out;
if ($outfile) {
open(my $fh, '>', $outfile) or die "Could not open file '$outfile' $!";
  $embl_out= Bio::SeqIO->new(-fh => $fh, -format => 'embl');
}
else{
  $embl_out = Bio::SeqIO->new(-fh => \*STDOUT, -format => 'embl');
}

### Read gb input file.
my $gb_in = Bio::SeqIO->new(-file => $gb, -format => 'genbank');


### MAIN ###

while( my $seq = $gb_in->next_seq) {

  $embl_out->write_seq($seq);

}

__END__

=head1 NAME

gaas_gb2embl.pl

=head1 DESCRIPTION

The script take a Genebank file as input, and will translate it in EMBL format.

=head1 SYNOPSIS

    gaas_gb2embl.pl --gb infile.gb [ -o outfile ]

=head1 OPTIONS

=over 8

=item B<--gb>

Input genebank file that will be read

=item B<-o> , B<--output> , B<--out> , B<--outfile> or B<--embl>

Output embl file. If no output file is specified, the output will be
written to STDOUT.

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
