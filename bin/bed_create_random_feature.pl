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

# PARAMETERS - OPTION
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
    pod2usage( { -verbose => 2,
                 -exitval => 0 } );
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

gffRepeat_analyzer.pl -
The script aims to create a fake bed file.

=head1 SYNOPSIS

    script.pl -g name -s 10000 -o <output file>
    script.pl --help

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

=cut
