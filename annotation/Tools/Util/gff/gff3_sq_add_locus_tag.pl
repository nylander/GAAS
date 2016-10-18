#!/usr/bin/env perl

###################################################
# Jacques Dainat 01/2016                          #  
# Bioinformatics Infrastructure for Life Sciences #
# jacques.dainat@bils.se                          #
###################################################

use Carp;
use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use IO::File ;
use Bio::Tools::GFF;
use BILS::Handler::GXFhandler qw(:Ok);
use BILS::Handler::GFF3handler qw(:Ok);

my $start_run = time();

my $inputFile=undef;
my $outfile=undef;
my $opt_help = 0;
my $locus_tag=undef;

Getopt::Long::Configure ('bundling');
if ( !GetOptions ('file|input|gff=s' => \$inputFile,
      'l=s' => \$locus_tag,
      'o|output=s' => \$outfile,
      'h|help!'         => \$opt_help )  )
{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

if ($opt_help) {
    pod2usage( { -verbose => 2,
                 -exitval => 0 } );
}

if ((!defined($inputFile)) ){
   pod2usage( { -message => 'at least 1 parameter is mandatory: -i',
                 -verbose => 1,
                 -exitval => 1 } );
}

my $ostream     = IO::File->new();

# Manage input fasta file
my $ref_in = Bio::Tools::GFF->new(-file => $inputFile, -gff_version => 3);

# Manage Output
my $gffout;
if ($outfile) {
  $outfile=~ s/.gff//g;
  open(my $fh, '>', $outfile.".gff") or die "Could not open file '$outfile' $!";
  $gffout= Bio::Tools::GFF->new(-fh => $fh, -gff_version => 3 );
  
}
else{
  $gffout = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 3);
}

#define the locus tag
if(! $locus_tag){
  $locus_tag="locus_tag";
}

#time to calcul progression
my $startP=time;
my $nbLine=`wc -l < $inputFile`;
$nbLine =~ s/ //g;
chomp $nbLine;
print "$nbLine line to process...\n";

my $geneName=undef;
my $line_cpt=0;
while (my $feature = $ref_in->next_feature() ) {
  $line_cpt++;

    if( lc($feature->primary_tag) eq "gene"){
      $geneName=$feature->_tag_value('ID');
      $feature->add_tag_value($locus_tag, $geneName);
    }
    elsif($geneName){
      $feature->add_tag_value($locus_tag, $geneName);
    }
    $gffout->write_feature($feature);

  #Display progression
  if ((30 - (time - $startP)) < 0) {
    my $done = ($line_cpt*100)/$nbLine;
    $done = sprintf ('%.0f', $done);
        print "\rProgression : $done % processed.\n";
    $startP= time;
  }
}

##Last round
my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job done in $run_time seconds\n";



__END__

=head1 NAME

gff3_add_locus_tag.pl -
add a locus tag based on the ID of the gene feature. 

=head1 SYNOPSIS

    gff3_add_locus_tag.pl --gff <input file> [-o <output file>]
    gff3_add_locus_tag.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>, B<--file> or B<--input>

STRING: Input gff file that will be read.

=item B<-l> 

Locus tag, by defaut it will be called locus_tag, but using this option you can specied the name of this attribute.

=item B<-o> or B<--output> 

STRING: Output file.  If no output file is specified, the output will be written to STDOUT. The result is in tabulate format.

=item B<--help> or B<-h>

Display this helpful text.

=back

=cut
