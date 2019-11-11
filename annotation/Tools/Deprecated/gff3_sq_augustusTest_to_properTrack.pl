#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use IO::File;
use List::MoreUtils qw(uniq);
use File::Basename;
use Bio::Tools::GFF;

my $header = qq{
########################################################
# NBIS 2014 - Sweden                                   #
# jacques.dainat\@nbis.se                               #
# Please cite NBIS (www.nbis.se) when using this tool. #
########################################################
};

my $gff = undef;
my $opt_output = undef;

my $help= 0;

if ( !GetOptions(
    "help|h" => \$help,
    'o|output=s'      => \$opt_output,
    "gff|f=s" => \$gff))

{
    pod2usage( { -message => "Failed to parse command line",
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($help) {
    pod2usage( { -verbose => 1,
                 -exitval => 0,
                 -message => "$header \n" } );
}

if ( ! (defined($gff)) ){
    pod2usage( {
           -message => "$header\nAt least 1 parameter is mandatory:\nInput reference gff file (--gff) \n\n",
           -verbose => 0,
           -exitval => 2 } );
}

#### OUT
my $gffout;
if ($opt_output) {
  ##### Stream out
  open($gffout, '>', $opt_output) or die "Could not open file '$opt_output' $!";
  }
else{
  open($gffout, '>&', STDOUT,) or die "Could not open file '$opt_output' $!";
}

                #####################
                #     MAIN          #
                #####################



######################
### Parse GFF input #
print "Reading file $gff\n";
my $fh1;
if ($gff) {
  open($fh1, '<', $gff) or die "Could not open file '$gff' $!";
}
### END Parse GFF input #
#########################

while( my $line = <$fh1>)  {
  my @list = split(/\s/,$line);
  my $header = $list[0];
  my @header_parts = split(/_/,$header);
  my @positions = split(/-/,$header_parts[1]);
  my $start = $positions[0];
  if($list[1] ne "database"){
    print $gffout $header_parts[0]."\t".$list[1]."\t".$list[2]."\t".($start+$list[3])."\t".($start+$list[4])."\t".$list[5]."\t".$list[6]."\t".$list[7]."\t".$list[8]."\n";
  }
}



print "Done\n";




__END__

=head1 NAME

gff3_sq_augustusTest_to_properTrack.pl -
The script take a gff3 file from Augustus  as input. It has to have been preproceed by gxf_to_gff.pl first.
It will recreate proper coordinate to visualise the gff file into a browser.

=head1 SYNOPSIS

    ./gff3_sq_augustusTest_to_properTrack.pl -gff file.gff  [ -o outfile ]
    ./gff3_sq_augustusTest_to_properTrack.pl --help

=head1 OPTIONS

=over 8

=item B<--gff> or B<-f>

Input GFF3 file that will be read (and sorted)

=item B<--output> or B<-o>

File where will be written the result. If no output file is specified, the output will be written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.

=back

=cut
