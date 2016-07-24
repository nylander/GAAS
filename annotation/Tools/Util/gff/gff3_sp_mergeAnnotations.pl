#!/usr/bin/perl


use Carp;
use Clone 'clone';
use strict;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use List::MoreUtils qw(uniq);
use Bio::Tools::GFF;
use BILS::Handler::GXFhandler qw(:Ok);
use BILS::Handler::GFF3handler qw(:Ok);

my $usage = qq{
########################################################
# BILS 2015 - Sweden                                   #	
# jacques.dainat\@bils.se                               #
# Please cite BILS (www.bils.se) when using this tool. #
########################################################

Usage: perl my_script.pl --gff1 Infile1 --gff2 Infile2 [--out outfile]
  Getting help:
    [--help]

  Input:
    [--gff1 filename1]
		The name of the gff3 file 1 used as reference.

    [--gff2 filename2]
    The name of the gff3 file 2. This file will be used to complete the file1.
	
  Ouput:    
    [--out filename]
        The name of the output file (A GFF file).

  This script merge two gff3 annotation together. If features exist in both file, only one is writen in output.

};

my $outfile = undef;
my $file1 = undef;
my $file2 = undef;
my $help= 0;

if ( !GetOptions(
    "help|h" => \$help,
    "gff1|f1=s" => \$file1,
    "gff2|f2=s" => \$file2,
    "output|outfile|out|o=s" => \$outfile))

{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($help) {
    pod2usage( { -verbose => 2,
                 -exitval => 0,
                 -message => "$usage\n" } );
}
 
if ( ! ((defined($file1)) and (defined($file2)))){
    pod2usage( {
           -message => "\nAt least 2 parameter is mandatory:\nInput reference gff file1 (--gff1) and Input reference gff file2 (--gff2)\n\n".
           "$usage\n",
           -verbose => 0,
           -exitval => 2 } );
}

######################
# Manage output file #
my $gffout;
if ($outfile) {
  $outfile=~ s/.gff//g;
open(my $fh, '>', $outfile.".gff") or die "Could not open file '$outfile' $!";
  $gffout= Bio::Tools::GFF->new(-fh => $fh, -gff_version => 3 );
}
else{
  $gffout = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 3);
}


                #####################
                #     MAIN          #
                #####################


######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) = BILS::Handler::GFF3handler->slurp_gff3_file_JD($file1);
print ("$file1 GFF3 file parsed\n");
info_omniscient($hash_omniscient);
my ($hash_omniscient2, $hash_mRNAGeneLink2) = BILS::Handler::GFF3handler->slurp_gff3_file_JD($file2);
print ("$file2 GFF3 file parsed\n");
info_omniscient($hash_omniscient2);
merge_omniscients($hash_omniscient2, $hash_omniscient);
print ("$file1 and $file2 merged:\n");
info_omniscient($hash_omniscient);

########
# Print results
print_omniscient($hash_omniscient, $gffout);  

__END__
