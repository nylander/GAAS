#!/usr/bin/env perl

## BILS 2015
## jacques.dainat@bils.se

use strict;
use Data::Dumper;
use Pod::Usage;
use Getopt::Long;
use Bio::SeqIO;

my $usage = qq{
########################################################
# BILS 2015 - Sweden                                   #  
# jacques.dainat\@bils.se                               #
# Please cite BILS (www.bils.se) when using this tool. #
########################################################

Usage: perl my_script.pl --gb Infile [--out outfile]
  Getting help:
    [--help]

  Input:
    [--gb filename]
    The name of the Genebank file to convert. 
  
  Ouput:    
    [--out filename]
        The name of the output file (An EMBL file).
};

my $outfile = undef;
my $gb = undef;
my $help;

if( !GetOptions(
    "help" => \$help,
    "gb=s" => \$gb,
    "outfile|output|o|out|embl=s" => \$outfile))
{
    pod2usage( { -message => "Failed to parse command line\n$usage",
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if ( ! (defined($gb)) ){
    pod2usage( {
           -message => "Missing the --gb argument\n$usage",
           -verbose => 0,
           -exitval => 2 } );
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
  #print $seq->primary_seq;
  #print Dumper($seq );
  $embl_out->write_seq($seq);
  exit;
}

__END__

=head1 NAME

gb2embl.pl -
The script take a Genebank file as input, and will translate it in EMBL format.

=head1 SYNOPSIS

    ./gb2embl.pl --gb=infile.gff [ -o outfile ]

=head1 OPTIONS

=over 8

=item B<-gb>

Input genebank file that will be read 

=item B<-o> , B<--output> , B<--out> , B<--outfile> or B<--embl>

Output embl file. If no output file is specified, the output will be
written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.

=back

=cut