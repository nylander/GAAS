#!/usr/bin/env perl

#############################################
# Jacques Dainat 2018
#############################################


#libraries
use strict;
use warnings;
use Getopt::Long;
use IO::File;
use Pod::Usage;
use Carp;
use Bio::FeatureIO;
use GAAS::GAAS;

# PARAMETERS - OPTION
my $header = get_gaas_header();
my $opt_genome;
my $opt_sizeGenome;
my $opt_sizeGeneMAx;
my $opt_nb;
my $opt_output;
my $opt_help;


# OPTION MANAGMENT
if ( !GetOptions( 'g|genome|fa=s' => \$opt_genome,
                  's|size=i' => \$opt_sizeGenome,
                  'nbg|number_gene=i' => \$opt_nb,
                  'sg|size_gene=i' => \$opt_sizeGeneMAx,
                  'o|output=s'      => \$opt_output,
                  'h|help!'         => \$opt_help ) )
{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

if ($opt_help) {
    pod2usage( { -verbose => 99,
                 -exitval => 0,
								 -message => "$header\n" } );
}

if ( ! (defined($opt_genome)) and ! (defined($opt_sizeGenome)) ){
    pod2usage( {
           -message => "\nAt least 1 parameter is mandatory:\nInput reference gff file (--f)\n\n".
           "Many optional parameters are available. Look at the help documentation to know more.\n",
           -verbose => 0,
           -exitval => 2 } );
}




if ( ! (defined($opt_sizeGeneMAx))) {
    print "you didnt define size of gene. We will use 1000 bp by default.\n";
    $opt_sizeGeneMAx=1000;
}
if ( ! (defined($opt_nb)) ){
    print "you didnt define number of gene. We will use 100 by default.\n";
    $opt_nb=100;
}

my $output;
if ($opt_output) {
      if (-f $opt_output){
      print "Cannot create a file with the name $opt_output because a file with this name already exists.\n";exit();
  }
    open(my $fh, '>', $opt_output) or die "Could not open file '$opt_output' $!";
     $output= Bio::FeatureIO->new(-fh => $fh, -format => 'BED' );
}
else{
  $output = Bio::FeatureIO->new(-fh => \*STDOUT, -format => 'BED' );
}


my $seq_id;
if ( ! (defined($opt_genome))) {
    $seq_id="chr_unknown";
}
else{
    $seq_id=$opt_genome;
}

for (my $i=0; $i <= $opt_nb; $i++) {
    my $start=int(rand($opt_sizeGenome-$opt_sizeGeneMAx));
    my $end=$start+$opt_sizeGeneMAx;

    my $primary_tag="gene_invent".$i;

    my $random_strand = int(rand(2));

    my $feature = Bio::SeqFeature::Annotated->new(-seq_id => $seq_id, -start => $start, -end => $end, -strand => $random_strand ) ;
    $output->write_feature($feature);
}
print "FINISH !!\n";


__END__

=head1 NAME

gaas_create_random_feature.pl

=head1 DESCRIPTION

The script aims to create a fake bed file.

=head1 SYNOPSIS

    gaas_create_random_feature.pl -g name -s 10000 -o <output file>
    gaas_create_random_feature.pl --help

=head1 OPTIONS

=over 8

=item B<-g>, B<--genome> or B<--fa>

STRING: Name to use for writing in first column of the bed file. default chr_unknown.

=item B<-s>, B<--size>

INTEGER: Genome size. It define the range where features will be created.

=item B<--nbg>, B<--number_gene>

INTEGER: Number of gene. It define the number of gene features to be created.

=item B<--sg>, B<--size_gene>

INTEGER: Size of genes. It define the size oft the gene features to be created.

=item B<-o> or B<--output>

STRING: Output file.  If no output file is specified, the output will be written to STDOUT. The result is in tabulate format.

=item B<--help> or B<-h>

Display this helpful text.

=back

=head1 FEEDBACK

=head2 Did you find a bug?

Do not hesitate to report bugs to help us keep track of the bugs and their
resolution. Please use the GitHub issue tracking system available at this
address:

            https://github.com/NBISweden/GAAS/issues

 Ensure that the bug was not already reported by searching under Issues.
 If you're unable to find an (open) issue addressing the problem, open a new one.
 Try as much as possible to include in the issue when relevant:
 - a clear description,
 - as much relevant information as possible,
 - the command used,
 - a data sample,
 - an explanation of the expected behaviour that is not occurring.

=head2 Do you want to contribute?

You are very welcome, visit this address for the Contributing guidelines:
https://github.com/NBISweden/GAAS/blob/master/CONTRIBUTING.md

=cut

AUTHOR - Jacques Dainat
