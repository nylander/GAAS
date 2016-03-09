#!/usr/bin/perl

## BILS 2015
## jacques.dainat@bils.se

use strict;
use Pod::Usage;
use Getopt::Long;
use POSIX qw(strftime);
use Bio::SeqIO;
use Data::Dumper;
use Bio::Tools::GFF;
use BILS::Handler::GXFhandler qw(:Ok);
use BILS::Handler::GFF3handler qw(:Ok);
use Bio::DB::Fasta;

my $usage = qq{
########################################################
# BILS 2015 - Sweden                                   #  
# jacques.dainat\@bils.se                               #
# Please cite BILS (www.bils.se) when using this tool. #
########################################################

Usage: perl my_script.pl --gff Infile [--out outfile]
  Getting help:
    [--help]

  Input:
    [--gff filename]
    The name of the GFF file to convert. 
  
    [--fasta filename]
    fasta file name.  

  Ouput:    
    [--out filename]
        The name of the output file (A EMBL file).
};

my $outfile = undef;
my $gff = undef;
my $file_fasta=undef;
my $help;

if( !GetOptions(
    "help" => \$help,
    "gff|in=s" => \$gff,
    "fasta|fa=s" => \$file_fasta,
    "outfile|output|o|out=s" => \$outfile))
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

if ( ! (defined($gff)) or ! (defined($file_fasta)) ){
    pod2usage( {
           -message => "Missing the --gff argument\n$usage",
           -verbose => 0,
           -exitval => 2 } );
}

## Manage output file
my $embl_out;
if ($outfile) {
open(my $fh, '>', $outfile) or die "Could not open file '$outfile' $!";
  $embl_out = Bio::SeqIO->new(-fh => $fh, -format => 'embl');
}
else{
  $embl_out = Bio::SeqIO->new(-fh => \*STDOUT, -format => 'embl');
}

######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) = BILS::Handler::GFF3handler->slurp_gff3_file_JD($gff);
print ("GFF3 file parsed\n");

my $hash_by_group=group_features_from_omniscient($hash_omniscient);
print ("GFF3 data grouped\n");
### MAIN ###

###
# or for GMT formatted appropriately for your locale:
my $datestring = strftime "%e-%h-%g", gmtime;

####################
# index the genome #
my $db = Bio::DB::Fasta->new($file_fasta);
print ("Genome fasta parsed\n");


foreach my $seq_id (keys %{$hash_by_group} ){
 
  my $seqObject = Bio::Seq->new(-seq => $db->seq($seq_id));
 # my $location = Bio::Location::Split->new();

  foreach my $gene_grouped(keys %{$hash_by_group->{$seq_id}}){


    foreach my $feature (@{$hash_by_group->{$seq_id}{$gene_grouped}}){
 #     $location->add_sub_Location($feature->location);
#      print $feature->location;
#      my $genef = Bio::SeqFeature::Generic->new(-location =>$location, -primary_tag => 'CDS');
      $seqObject->add_SeqFeature($feature);
    } 
  }
  print Dumper($seqObject);
  $embl_out->write_seq($seqObject);
}

#my $db->seq( $seqid)
#while( my $feature = $gtfio->next_feature) {
#  $embl_out->write_seq($seq)
#}

__END__

It is probably easiest to just group things and make a split location.  
You will have the most control over the objects you create.

my %genes;
while( my $f = $gff->next_feature ) {
   my ($group) = $feature->get_tag_values('Group'); # substitute group 
with whatever you have in the group field
  push @{$gene{$group}}, $feature;
}
# get a Bio::Seq object called $seq somehow, either by reading in a 
fasta sequence file, etc...
while( my ($gene,$features) = each %genes ) {
  my $location = Bio::Location::Split->new();
  for my $f ( @$features ) {
    $location->add_sub_Location($f->location);
  }
  my $genef = Bio::SeqFeature::Generic->new(-location =>$location, 
-primary_tag => 'CDS');
  $seq->add_SeqFeature($genef);
}
my $seqio = Bio::SeqIO->new(-format => 'genbank');
$seqio->write_seq($seq);

=head1 NAME

gb2embl.pl -
The script take a EMBL file as input, and will translate it in Genbank format.

=head1 SYNOPSIS

    ./embl2gb.pl --embl=infile.gff [ -o outfile ]

=head1 OPTIONS

=over 8

=item B<--embl>

Input EMBL file that will be read 

=item B<-o> , B<--output> , B<--out> , B<--outfile> or B<--gff>

Output GFF file. If no output file is specified, the output will be
written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.

=back

=cut