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
use Bio::DB::Fasta;
use IO::File ;
use Bio::Tools::GFF;
use BILS::Handler::GXFhandler qw(:Ok);
use BILS::Handler::GFF3handler qw(:Ok);

my $start_run = time();

my $opt_gfffile=undef;
my $opt_fastafile=undef;
my $outfile=undef;
my $opt_help = 0;


Getopt::Long::Configure ('bundling');
if ( !GetOptions ('file|input|gff=s' => \$opt_gfffile,
      'f|fasta=s' => \$opt_fastafile,
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

if ((!defined($opt_gfffile)) ){
   pod2usage( { -message => 'at least 2 parameters are mandatory',
                 -verbose => 1,
                 -exitval => 1 } );
}

my $ostream     = IO::File->new();

# Manage input fasta file
my $ref_in = Bio::Tools::GFF->new(-file => $opt_gfffile, -gff_version => 3);

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


#### read fasta
my $nbFastaSeq=0;
my $db = Bio::DB::Fasta->new($opt_fastafile);
print ("Genome fasta parsed\n");

#time to calcul progression
my $startP=time;
my $cpt_removed=0;
while (my $feature = $ref_in->next_feature() ) {

  if($db->seq($feature->seq_id)){
    #create sequence object
    $feature->gff_string($gffout) 
  }
  else{
    print "SequenceID ".$feature->seq_id." is absent from the fasta file\n";
    $cpt_removed++;
  }
}

print "We removed $cpt_removed annotation.\n";
my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job done in $run_time seconds\n";



__END__

=head1 NAME

gff3_sq_keep_annotation_from_fastaSeq.pl -
This script is a kind of annotation filter by sequence name. It goes through the gff annotation features and remove those that are not linked to a sequence from the fasta file provided.

=head1 SYNOPSIS

    gff3_sq_keep_annotation_from_fastaSeq.pl --gff <gff_file.gff> --fasta <fasta_file.fa> [-o <output file>]
    gff3_sq_keep_annotation_from_fastaSeq.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>, B<--file> or B<--input>

STRING: Input gff file.

=item B<-f> or B<--fasta>

STRING: fasta file.

=item B<-o> or B<--output> 

STRING: Output file.  If no output file is specified, the output will be written to STDOUT. The result is in tabulate format.

=item B<--help> or B<-h>

Display this helpful text.

=back

=cut
