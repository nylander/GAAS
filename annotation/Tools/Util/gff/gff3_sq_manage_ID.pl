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
use BILS::Handler::GFF3handler qw(:Ok);
use BILS::Handler::GXFhandler qw(:Ok);

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
my %hash_IDs;
my %featCount;
my %mapID;
while (my $feature = $ref_in->next_feature() ) {
  $line_cpt++;

  _uniq_ID ($feature, \%hash_IDs, \%featCount, \%mapID);
  
  if($feature->has_tag('Parent')){
    my $parent = lc($feature->_tag_value('Parent'));
    if(! exists($mapID{$parent})){
      print "How is it possible ? This parent hasn't been seen before\n";
    }
     create_or_replace_tag($feature,'Parent', $mapID{$parent});
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



sub _uniq_ID{
  my ($feature, $hash_IDs, $miscCount, $mapID) = @_;

  
  my  $key=lc($feature->primary_tag);
  $miscCount->{$key}++;
  my $id = $key."-".$miscCount->{$key};

  while( exists_keys($hash_IDs, ($id) ) ){  #loop until we found an uniq tag 
    $miscCount->{$key}++;
    $id = $key."-".$miscCount->{$key};
  }

  #push the new ID  
  $hash_IDs->{$id}++;
  $mapID->{lc($feature->_tag_value('ID'))} = $id;

  # modify the feature ID with the correct one chosen
  create_or_replace_tag($feature,'ID', $id); #modify ID to replace by parent value
}

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
