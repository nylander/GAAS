#!/usr/bin/env perl

###################################################
# Jacques Dainat 01/2016                          #
# Bioinformatics Infrastructure for Life Sciences #
# jacques.dainat@nbis.se                          #
###################################################

use strict;
use warnings;
use Carp;
use Pod::Usage;
use Getopt::Long;
use IO::File ;
use Bio::Tools::GFF;
use NBIS::GFF3::Omniscient qw(select_gff_format create_or_replace_tag);

my $start_run = time();

my $inputFile=undef;
my $outfile=undef;
my $outformat=undef;
my $opt_help = 0;

Getopt::Long::Configure ('bundling');
if ( !GetOptions ('file|input|gff|i=s' => \$inputFile,
      'of=i' => \$outformat,
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


# Manage input fasta file
my $format = select_gff_format($inputFile);
my $ref_in = Bio::Tools::GFF->new(-file => $inputFile, -gff_version => $format);


# Manage Output
if(! $outformat){
  $outformat=$format;
}

my $gffout;
if ($outfile) {
  $outfile=~ s/.gff//g;
  open(my $fh, '>', $outfile.".gff") or die "Could not open file '$outfile' $!";
  $gffout= Bio::Tools::GFF->new(-fh => $fh, -gff_version => $outformat );

}
else{
  $gffout = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => $outformat);
}

#time to calcul progression
my $startP=time;
my $nbLine=`wc -l < $inputFile`;
$nbLine =~ s/ //g;
chomp $nbLine;
print "$nbLine line to process...\n";

my $line_cpt=0;
my $previousGeneID="";
my $gene_id=0;
my $previousTranscriptID="";
my $transcript_id=0;
while (my $feature = $ref_in->next_feature() ) {
  $line_cpt++;

  if($feature->has_tag('gene_id')){
    my $parent = lc($feature->_tag_value('gene_id'));
    if($parent ne $previousGeneID){
      $gene_id++;
      $previousGeneID=$parent;
    }
     create_or_replace_tag($feature,'gene_id', $gene_id);
  }
  if($feature->has_tag('transcript_id')){
    my $parent = lc($feature->_tag_value('transcript_id'));
    if($parent ne $previousTranscriptID){
      $transcript_id++;
      $previousTranscriptID=$parent;
    }
     create_or_replace_tag($feature,'transcript_id', $gene_id);
  }

  $gffout->write_feature($feature);

  #####################
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

gff3_sq_manage_ID.pl -
change IDs to give uniq one. This script is sequential, it means it will works great even if 2 different group of feature an parent that have the same ID. At the end they will have different IDs.

=head1 SYNOPSIS

    gff3_sq_manage_ID.pl --gff <input file> [-o <output file>]
    gff3_sq_manage_ID.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>, B<--file> or B<--input>

STRING: Input gff(1 or 2 or 3) or gtf file that will be read.

=item B<--of>

Output format, if no ouput format is given, the same as the input one detected will be used. Otherwise you can force to have a gff version 1 or 2 or 3 by giving the corresponding number.

=item B<-o> or B<--output>

STRING: Output file.  If no output file is specified, the output will be written to STDOUT. The result is in tabulate format.

=item B<--help> or B<-h>

Display this helpful text.

=back

=cut
