#!/usr/bin/perl

####
# Jacques Dainat 2015/03
# jacques.dainat@bils.se
####
use strict;
use Pod::Usage;
use Getopt::Long;
use Bio::SeqIO ;
use Bio::DB::Fasta;
use Bio::Tools::GFF;
use BILS::Handler::GXFhandler qw(:Ok);
use BILS::Handler::GFF3handler qw(:Ok);


my $start_run = time();
my $opt_gfffile;
my $opt_fastafile;
my $opt_output;
my $opt_AA=undef;
my $opt_help = 0;

my $width = 60; # line length printed

my $header = qq{
########################################################
# BILS 2015 - Sweden                                   #  
# jacques.dainat\@bils.se                               #
# Please cite BILS (www.bils.se) when using this tool. #
########################################################
};

# OPTION MANAGMENT
if ( !GetOptions( 'g|gff=s' => \$opt_gfffile,
                  'f|fa|fasta=s' => \$opt_fastafile,
                  'protein|p|aa' => \$opt_AA,
                  'o|output=s'      => \$opt_output,
                  'h|help!'         => \$opt_help ) )
{
    pod2usage( { -message => "$header\nFailed to parse command line",
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($opt_help) {
    pod2usage( { -verbose => 2,
                 -exitval => 2,
                 -message => "$header \n" } );
}
 
if ( (! (defined($opt_gfffile)) ) || (! (defined($opt_fastafile)) ) ){
    pod2usage( {
           -message => "\nAt least 2 parametes are mandatory:\nInput reference gff file (-g);  Input reference fasta file (-f)\n\n".
           "Output is optional. Look at the help documentation to know more.\n",
           -verbose => 0,
           -exitval => 2 } );
}


my $ostream;
if ($opt_output) {
  $opt_output=~ s/.fasta//g;
  $opt_output=~ s/.fa//g;
  open(my $fh, '>', $opt_output.".fa") or die "Could not open file '$opt_output' $!";
  $ostream= Bio::SeqIO->new(-fh => $fh, -format => 'Fasta' );
}
else{
  $ostream = Bio::SeqIO->new(-fh => \*STDOUT, -format => 'Fasta');
}

##### MAIN ####

#### read gff file and save info in memory
######################
### Parse GFF input #
print "Reading file $opt_gfffile\n";
my ($hash_omniscient, $hash_mRNAGeneLink) = BILS::Handler::GFF3handler->slurp_gff3_file_JD($opt_gfffile);
print "Parsing Finished\n";
### END Parse GFF input #
#########################

my $hash_l1_grouped = group_l1features_from_omniscient($hash_omniscient);

#### read fasta
my $nbFastaSeq=0;
#my $inFasta  = Bio::SeqIO->new(-file => "$opt_fastafile" , '-format' => 'Fasta');
my $db = Bio::DB::Fasta->new($opt_fastafile);
print ("Genome fasta parsed\n");

  foreach my $seqname (keys $hash_l1_grouped) {
     foreach my $feature_l1 (@{$hash_l1_grouped->{$seqname}}) {
        my $id_l1=lc($feature_l1->_tag_value('ID'));
        my $name=undef;

        if ($feature_l1->has_tag('Name')){
          $name = $feature_l1->_tag_value('Name');
        }
        elsif($feature_l1->has_tag('gene')){
          $name = $feature_l1->_tag_value('gene');
        }


        #################
        # == LEVEL 2 == #
        #################
        foreach my $ptag_l2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
          
          if ( exists ($hash_omniscient->{'level2'}{$ptag_l2}{$id_l1} ) ){
            foreach my $feature_l2 ( @{$hash_omniscient->{'level2'}{$ptag_l2}{$id_l1}}) {

              my $id_l2  = lc($feature_l2->_tag_value('ID'));

              my $header=$id_l2;
              if ($feature_l2->has_tag('Name') and ! $name){
                $name = $feature_l2->_tag_value('Name');
              }
              elsif($feature_l2->has_tag('gene') and ! $name){
                $name = $feature_l2->_tag_value('gene');
              }

              #################
              # == LEVEL 3 == #
              #################

              if ( exists ($hash_omniscient->{'level3'}{'cds'}{$id_l2} ) ){
                my $seqObj = extract_cds_sequence($hash_omniscient->{'level3'}{'cds'}{$id_l2}, $db);

                my $header=$id_l2;
                if($name){
                  $header.="|Name=".$name;
                }
                $seqname =~ tr/|/_/;
                $header.="|Seq_id=".$seqname;
                $seqObj->id($header);

                if($opt_AA){ #translate if asked
                  my $transObj = $seqObj->translate();
                  
                  $ostream->write_seq($transObj);  
                }
                else{
                  $ostream->write_seq($seqObj);                
                }
                $nbFastaSeq++;
              }
            }
          }
        }
      }
    }

print "$nbFastaSeq cds converted in fasta.\n";
my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job done in $run_time seconds\n";

__END__

=head1 NAME

gff3_extract_cds_sequences.pl -
This script ... 
The result is written to the specified output file, or to STDOUT.

=head1 SYNOPSIS

    ./gff3_extract_cds_sequences.pl -g=infile.gff -f=infile.fasta  [ -o outfile ]
    ./gff3_extract_cds_sequences.pl --help

=head1 OPTIONS

=over 8

=item B<-g>, B<--gff> or B<-ref>

Input GFF3 file that will be read (and sorted)

=item B<-f> or B<--fasta> 

Input fasta file that will be masked

=item B<-sm> 

SoftMask option =>Sequences masked will be in lowercase

=item B<-hm> 

HardMask option => Sequences masked will be replaced by a character. By default the character used is 'n'. But you are allowed to speceify any character of your choice. To use 'z' instead of 'n' type: -hm z

=item B<-o> or B<--output>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.

=back

=cut
