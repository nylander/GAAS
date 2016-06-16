#!/usr/bin/perl

###################################################
# domainExtractor.pl - Jacques Dainat 01/2015     #
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

my $inputFile;
my $outputFile;
my $genomeSize;
my $opt_help = 0;

Getopt::Long::Configure ('bundling');
if ( !GetOptions ('i|file|input|gff=s' => \$inputFile,
      'o|output=s' => \$outputFile,
      'g|genome=i' => \$genomeSize,
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
if(defined($outputFile))
{
$ostream->open( $outputFile, 'w' ) or
  croak(
     sprintf( "Can not open '%s' for reading: %s", $outputFile, $! ) );
}
else{
  $ostream->fdopen( fileno(STDOUT), 'w' ) or
      croak( sprintf( "Can not open STDOUT for writing: %s", $! ) );
}

#time to calcul progression
my $startP=time;print "start\n";
my $nbLine=`wc -l $inputFile | cut -d' ' -f4`;
chomp $nbLine;
print "$nbLine line to process...\n";
my $line_cpt=0;

my $type_count;
my $type_bp;
while (my $feature = $ref_in->next_feature() ) {
  $line_cpt++;
  my $type = lc($feature->primary_tag);
  ## repeatMasker or repeatRunner
  if (($type eq 'match') or ($type eq 'protein_match')){

     my $nameAtt=$feature->_tag_value('Name');
     my $genus=(split ":", (split /\|/, (split /\s+/,$nameAtt)[0])[-1])[-1];
     $type_count->{$genus}++;
     $type_bp->{$genus}+=($feature->end()-$feature->start())+1;
  }

  #Display progression
  if ((1 - (time - $startP)) < 0) {
    my $done = ($line_cpt*100)/$nbLine;
    $done = sprintf ('%.0f', $done);
        print "Progression : $done % processed.\n";
    $startP= time;
  }
}
my $totalNumber=0;
my $totalSize=0;

if(defined($genomeSize)){
print $ostream "Repeat type\tNumber\tSize total (kb)\tSize mean (bp)\t% of the genome\t/!\\Results are rounding to two decimal places \n";
  foreach my $gnx (keys(%$type_count)) {
    my $Sitotal=sprintf("%0.2f",($type_bp->{$gnx}/1000));
    my $SizeMean=sprintf("%0.2f",($type_bp->{$gnx}/$type_count->{$gnx}));
    my $xGenome=sprintf("%0.2f",($type_bp->{$gnx}/$genomeSize)*100);
    print $ostream $gnx,"\t",$type_count->{$gnx},"\t",$Sitotal,"\t",$SizeMean,"\t",$xGenome,"\n";
    
    $totalNumber=$totalNumber+$type_count->{$gnx};
    $totalSize=$totalSize+$type_bp->{$gnx};

  }
}
else{
  print $ostream "Repeat type\tNumber\tSize total (kb)\tSize mean (bp)\t/!\\Results are rounding to two decimal places \n";
  foreach my $gnx (keys(%$type_count)) {
    my $Sitotal=sprintf("%0.2f",($type_bp->{$gnx}/1000));
    my $SizeMean=sprintf("%0.2f",($type_bp->{$gnx}/$type_count->{$gnx}));
    print $ostream $gnx,"\t",$type_count->{$gnx},"\t",$Sitotal,"\t",$SizeMean,"\n";

    $totalNumber=$totalNumber+$type_count->{$gnx};
    $totalSize=$totalSize+$type_bp->{$gnx};

  }
}

my $goodTotalSize=sprintf("%0.2f",($totalSize/1000));
my $goodTotalSizeMean=sprintf("%0.2f",($totalSize/$totalNumber));

if(defined($genomeSize)){
  my $goodxGenome=sprintf("%0.2f",($totalSize/$genomeSize)*100);
  print $ostream "Total\t",$totalNumber,"\t",$goodTotalSize,"\t",$goodTotalSizeMean,"\t",$goodxGenome,"\n";
}
else{
  print $ostream "Total\t",$totalNumber,"\t",$goodTotalSize,"\t",$goodTotalSizeMean,"\n";
}

my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job done in $run_time seconds\n";

__END__

=head1 NAME

gffRepeat_analyzor.pl -
The script allows to generate a tabulated format report of repeats annotated from a gff file containing repeats. 

=head1 SYNOPSIS

    gffRepeat_analyzor.pl -i <input file> [-g <integer> -o <output file>]
    gffRepeat_analyzor.pl --help

=head1 OPTIONS

=over 8

=item B<-i>, B<--gff>, B<--file> or B<--input>

STRING: Input gff file that will be read.

=item B<-g>, B<--genome>

INTEGER: Input genome size. It allows to calculate the percentage of the genome represented by each kind of repeats.

=item B<-o> or B<--output> 

STRING: Output file.  If no output file is specified, the output will be written to STDOUT. The result is in tabulate format.


=back

=cut
