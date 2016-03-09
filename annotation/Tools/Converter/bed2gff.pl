#!/usr/bin/env perl

## BILS 2015
## jacques.dainat@bils.se
# BED format described here: https://genome.ucsc.edu/FAQ/FAQformat.html#format1

use strict;
use Pod::Usage;
use Getopt::Long;
use Bio::Tools::GFF;
use Data::Dumper;

my $usage = qq{
########################################################
# BILS 2015 - Sweden                                   #  
# jacques.dainat\@bils.se                               #
# Please cite BILS (www.bils.se) when using this tool. #
########################################################
/!\\Only first 6 column inplemented... if your bed file contains more columns and you need their information... you need to finish the implementation };

my $outfile = undef;
my $bed = undef;
my $help;

if( !GetOptions(
    "help" => \$help,
    "bed=s" => \$bed,
    "outfile|output|o|out|gff=s" => \$outfile))
{
    pod2usage( { -message => "$header\nFailed to parse command line.\n",
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($help) {-message => "$header", 
            -verbose => 2,
            -exitval => 2 } );
}

if ( ! (defined($bed)) ){
    pod2usage( {
           -message => "$header\nMissing the --bed argument.\n",
           -verbose => 0,
           -exitval => 1 } );
}

## Manage output file
my $gffout;
if ($outfile) {
open(my $fh, '>', $outfile) or die "Could not open file '$outfile' $!";
  $gffout= Bio::Tools::GFF->new(-fh => $fh, -gff_version => 3);
}
else{
  $gffout = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 3);
}

# Ask for specific GFF information
print "Some information needed in a GFF3 file dont exist in a BED file... as we cannot guess is, please fill the following information:\n\n",
"Enter a source (3rd field in a gff file). Example: Cufflinks,Maker,Augustus,etc. [default: data]\n";
my $source_tag = <STDIN>;
chomp $source_tag;
if ($source_tag eq '') {$source_tag = 'data';}
if ($source_tag =~ /\s/) {die("Can't have whitespace in $source_tag\n") }

print "What is the data type ? Example: gene,mRNA,CDS,etc.  [default: match]\n";
my $primary_tag = <STDIN>;
chomp $primary_tag;
if ($primary_tag eq '') {$primary_tag = 'match';}
if ($primary_tag =~ /\s/) {die("Can't have whitespace in $primary_tag\n") }


### Read bed input file.
open my $fh, $bed or die "Could not open $bed: $!";






                                #######################
                                #        MAIN         #
                                #######################

my %bedOmniscent;
my $UniqID=0;
my $cpt_warning=0;
while( my $line = <$fh>)  {   
    chomp $line;

  if ($line =~ /#/){next;} #skip commented lines
    
  if (! $line =~ /\t/) {die("$line <> is not a tabulated format !\n") }
  else{

      my @fields = split /\t/, $line;
      if ($#fields == 0){
        if($cpt_warning == 0){
           print "This file doesnt look tabulated. BAD BOY ! So I will try to continue with space sparator.\n";
           $cpt_warning++;
        } 
        @fields = split /\s/, $line;
      }
      my $fieldNumber=$#fields+1;
      if($fieldNumber < 3 or $fieldNumber >12){
        print "Problem with that line:\n$line\nA bed file has at least three required fields ! 9 others fields are optional. So, a maximum of 12 fields is allowed !",
        "\n Your line contain $fieldNumber fields. Check the sanity of your file. BYE.\n";exit;
      }
      if($fieldNumber >= 7){
        print "sorry I have to implement to take in account all the optional fields !\n";exit;
      }
      my $cptField=0;
      $UniqID++;
      foreach my $field (@fields){
        $cptField++;

        ##########
        #MANDATORY fields
        ###########

        # 1 chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
        if($cptField == 1){
          $bedOmniscent{$UniqID}{'chrom'}=$field;
        }

        # 2 chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
        if($cptField == 2){
          $bedOmniscent{$UniqID}{'chromStart'}=$field;
        }

        # 3 chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
        if($cptField == 3){
          $bedOmniscent{$UniqID}{'chromEnd'}=$field;
        }
              
        ##########
        # OPTIONAL fields
        ##########

        # 4 name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
        if($cptField == 4){
          $bedOmniscent{$UniqID}{'name'}=$field;
        }
        
        # 5 score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). 
        if($cptField == 5){
          $bedOmniscent{$UniqID}{'score'}=$field;
        }

        # 6 strand - Defines the strand - either '+' or '-'.
        if($cptField == 6){
          $bedOmniscent{$UniqID}{'strand'}=$field;
        }
        
###==== I STOPPED HERE the implemetation (need to spend time on)

        # 7 thickStart - The starting position at which the feature is drawn thickly (for example, the start codon in gene displays). When there is no thick part, thickStart and thickEnd are usually set to the chromStart position.
        if($cptField == 7){
          $bedOmniscent{$UniqID}{'thickStart'}=$field;
        }

        # 8 thickEnd - The ending position at which the feature is drawn thickly (for example, the stop codon in gene displays).
        if($cptField == 8){
          $bedOmniscent{$UniqID}{'thickEnd'}=$field;
        }

        # 9 itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.              
        if($cptField == 9){
          $bedOmniscent{$UniqID}{'itemRgb'}=$field;
        }
        
        # 10 blockCount - The number of blocks (exons) in the BED line.
        if($cptField == 10){
          $bedOmniscent{$UniqID}{'blockCount'}=$field;
        }
        
        # 11 blockSizes - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
        if($cptField == 11){
          $bedOmniscent{$UniqID}{'blockSizes'}=$field;
        }
        
        # 12 blockStarts - A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount.
        if($cptField == 12){
          $bedOmniscent{$UniqID}{'blockStarts'}=$field;
        }
      }
  }
}


###
# MANAGE BED OMNISCIEN FOR OUTPUT
foreach my $id (keys %bedOmniscent){
#  foreach my $key (keys %{$bedOmniscent{$id}}){

    my $seq_id=$bedOmniscent{$id}{'chrom'};
    #my $source_tag; #fill at the beginning
    #my $primary_tag; #fill at the beginning
    my $start=$bedOmniscent{$id}{'chromStart'};
    my $end=$bedOmniscent{$id}{'chromEnd'};
    my $frame;

    my $score;
    if(exists($bedOmniscent{$UniqID}{'score'})){
      $score=$bedOmniscent{$UniqID}{'score'};
    }

    my $strand;
    if(exists($bedOmniscent{$UniqID}{'strand'})){
      $strand=$bedOmniscent{$UniqID}{'strand'};
    }

    my $feature = Bio::SeqFeature::Generic->new(-seq_id => $seq_id, -source_tag => $source_tag, -primary_tag => $primary_tag, -start => $start,  -end => $end , -frame => $frame , -strand =>$strand, tag => {'ID' => $id} ) ;
    
    if(exists($bedOmniscent{$id}{'name'})){
      $feature->add_tag_value('Name',$bedOmniscent{$id}{'name'});
    }

    $gffout->write_feature($feature);
}


close $fh;


__END__

=head1 NAME

bed2gff.pl -
The script take a bed file as input, and will translate it in gff format. /!\ Not implemented in it's totality...

=head1 SYNOPSIS

    ./bed2gff.pl --bed=infile.bed [ -o outfile ]

=head1 OPTIONS

=over 8

=item B<--bed>

Input bed file that will be convert.

=item B<-o> , B<--output> , B<--out> , B<--outfile> or B<--gff>

Output GFF file. If no output file is specified, the output will be
written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.

=back

=cut