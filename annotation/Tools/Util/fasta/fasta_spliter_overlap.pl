#!/usr/bin/env perl

###
# Implement case insensitive
###
use strict;
use Pod::Usage;
use Getopt::Long;
use Bio::SeqIO ;

my $start_run = time();

my $opt_fastafile;
my $opt_output;
my $opt_help = 0;
my $opt_chunck_size = undef;
my $opt_overlap = 0; 

my $header = qq{
########################################################
# BILS 2018 - Sweden                                   #  
# jacques.dainat\@nbis.se                               #
# Please cite NBIS (www.nbis.se) when using this tool. #
########################################################
};

# OPTION MANAGMENT
my @copyARGV=@ARGV;
if ( !GetOptions( 'f|fa|fasta=s' => \$opt_fastafile,
                  'c|chunck_size=s' => \$opt_chunck_size,
                  'l|overlap=s' => \$opt_overlap,
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
 
if (! defined($opt_fastafile) or ! defined($opt_chunck_size) ) {
    pod2usage( {
           -message => "\nAt least 2 parameter is mandatory:\nInput reference fasta file (-f)\nChunck_size (-c)\n".
           "Output is optional. Look at the help documentation to know more.\n",
           -verbose => 0,
           -exitval => 2 } );
}


#D OUTPUT
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

print "We will split fasta sequences by chunck of $opt_chunck_size and overlap of $opt_overlap.\n";
if($opt_overlap >=  $opt_chunck_size){
  print "$opt_overlap cannot be >= to $opt_chunck_size otherwise we end up in an infinite loop.\n";
}


##### MAIN ####

######### read fasta file #############
my $fasta1  = Bio::SeqIO->new(-file => $opt_fastafile , -format => 'Fasta');
while ( my $seq = $fasta1->next_seq() ) {
  my $start = 1;
  my $end = $opt_chunck_size;
  
  while ( $end < $seq->length() ) {

      my $sequence = undef;
      my $seqObj = undef;
      my $id_seq = undef;
    	if($seq->length() > ($end+$opt_chunck_size) ){
        
        
        $sequence = $seq->subseq($start, $end);
	$seqObj = Bio::Seq->new( '-format' => 'fasta' , -seq => $sequence);
        $id_seq = $seq->id."_".$start."_".$end;
        $seqObj->id($id_seq);
        $ostream->write_seq($seqObj);
    	}
      else{
        $sequence = $seq->subseq($start, $seq->length());
        $seqObj  = Bio::Seq->new( '-format' => 'fasta' , -seq => $sequence);
        $id_seq = $seq->id."_".$start."_".$seq->length();
        $seqObj->id($id_seq);
        $ostream->write_seq($seqObj);
	     last;
      }
      $start = $end - $opt_overlap;
      $end = $start + $opt_chunck_size;
  }
}


#END
print "usage: $0 @copyARGV\n";
my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job done in $run_time seconds\n";

#######################################################################################################################
        ####################
         #     METHODS    #
          ################
           ##############
            ############
             ##########
              ########
               ######
                ####
                 ##          




__END__

=head1 NAME

fasta_spliter_overlap.pl -
This script split sequences by size with an overlaped part.

=head1 SYNOPSIS

    ./fasta_filer_by_size.pl -f=infile.fasta [ -o outfile ]
    ./fasta_filer_by_size.pl --help

=head1 OPTIONS

=over 8

=item B<-f> or B<--fasta> 

Input fasta file.

=item B<-s>, B<--size>

Integer corresponding to a size in bp. Default value 1000. Sequence under the value will be discarded from the output.

=item B<-o> or B<--output>

Output fasta file.  If no output file is specified, the output will be
written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.

=back

=cut
