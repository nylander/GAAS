#!/usr/local/bin/perl -w

## BILS 2015
## jacques.dainat@bils.se

use strict;
use Pod::Usage;
use Getopt::Long;
use Bio::SeqIO;

my $usage = qq{
########################################################
# BILS 2015 - Sweden                                   #  
# jacques.dainat\@bils.se                               #
# Please cite BILS (www.bils.se) when using this tool. #
########################################################

Usage: perl my_script.pl --embl Infile [--out outfile]
  Getting help:
    [--help]

  Input:
    [--embl filename]
    The name of the EMBL file to convert. 
  
  Ouput:    
    [--out filename]
        The name of the output file (A Genbank file).
};

my $outfile = undef;
my $embl = undef;
my $help;

if( !GetOptions(
    "help" => \$help,
    "embl=s" => \$embl,
    "outfile|output|o|out|gb=s" => \$outfile))
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
           -message => "Missing the --embl argument\n$usage",
           -verbose => 0,
           -exitval => 2 } );
}

## Manage output file
my $gb_out;
if ($outfile) {
open(my $fh, '>', $outfile) or die "Could not open file '$outfile' $!";
  $gb_out= Bio::SeqIO->new(-fh => $fh, -format => 'genbank');
}
else{
  $gb_out = Bio::SeqIO->new(-fh => \*STDOUT, -format => 'genbank');
}

### Read gb input file. 
my $embl_in = Bio::SeqIO->new(-file => $embl, -format => 'embl');


### MAIN ###

while( my $seq = $embl_in->next_seq) {
  $gb_out->write_seq($seq)
}

__END__

=head1 NAME

gb2embl.pl -
The script take a EMBL file as input, and will translate it in Genbank format.

=head1 SYNOPSIS

    ./embl2gb.pl --embl=infile.gff [ -o outfile ]

=head1 OPTIONS

=over 8

=item B<-embl>

Input EMBL file that will be read 

=item B<-o> , B<--output> , B<--out> , B<--outfile> or B<--gb>

Output Genbank file. If no output file is specified, the output will be
written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.

=back

=cut