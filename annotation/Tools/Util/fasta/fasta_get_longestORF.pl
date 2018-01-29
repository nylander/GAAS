#!/usr/bin/env perl

use Carp;
use strict;
use POSIX qw(strftime);
use Getopt::Long;
use BILS::FASTA::Longest_orf;
use Pod::Usage;
use Bio::Seq;
use Bio::SeqIO;

my $start_run = time();

my $header = qq{
########################################################
# BILS 2017 - Sweden                                   #  
# jacques.dainat\@bils.se                               #
# Please cite BILS (www.bils.se) when using this tool. #
########################################################
};

my $codonTableId=1;
my $MIN_PROT_LENGTH = 100;
my $force_start_codon = undef;
my $force_complete = undef;
my $file_fasta=undef;
my $keep_all_orf=undef;
my $outfile = undef;
my $verbose = undef;
my $help= 0;

my @copyARGV=@ARGV;
Getopt::Long::Configure ('bundling');
if ( !GetOptions(
    "help|h" => \$help,
    "fasta|fa|f=s" => \$file_fasta,
    "size_min|s=i" => \$MIN_PROT_LENGTH,
    "force_start_codon!" => \$force_start_codon,
    "force_complete!" => \$force_complete,
    "table|codon|ct=i" => \$codonTableId,
    "keep_all_orf!" => \$keep_all_orf,
    "v!" => \$verbose,
    "output|outfile|out|o=s" => \$outfile))

{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($help) {
    pod2usage( { -verbose => 2,
                 -exitval => 2,
                 -message => "$header\n" } );
}
 
if ( !(defined($file_fasta)) ){
    pod2usage( {
           -message => "$header\nAt least 1 parameter is mandatory:\n Input fasta file (--fasta)\n\n",
           -verbose => 0,
           -exitval => 1 } );
}

if($codonTableId<0 and $codonTableId>25){
  print "$codonTableId codon table is not a correct value. It should be between 0 and 25 (0,23 and 25 can be problematic !)\n";
}

######################
# Manage output file #
my $fasta_out;
if ($outfile) {
  $fasta_out = Bio::SeqIO->new(-file => ">$outfile" , -format => 'fasta');
}
else{
  $fasta_out = Bio::SeqIO->new(-fh => \*STDOUT, -format => 'fasta');
}

# print usage performed
my $stringPrint = strftime "%m/%d/%Y at %Hh%Mm%Ss", localtime;
$stringPrint = "Launched the ".$stringPrint."\nusage: $0 @copyARGV\n";
print $stringPrint;

                #####################
                #     MAIN          #
                #####################

# .. Read genome fasta file.
my $inseq = Bio::SeqIO->new(-file   => "<$file_fasta", -format => 'fasta');

my %canditates;
while( my $seqObj = $inseq->next_seq() ) {


  my $longest_orf_finder = Longest_orf->new();
  $longest_orf_finder->allow_5prime_partials();
  $longest_orf_finder->allow_3prime_partials();

  my $seq = $seqObj->seq();
  my @orf_structs = $longest_orf_finder->capture_all_ORFs($seq);

  #sorting ORF by size
  @orf_structs = reverse sort {$a->{length}<=>$b->{length}} @orf_structs;
  my $acc = $seqObj->id();
  print "looking at sequence $acc \n" if $verbose;
  
  my %candidates;        
  while (@orf_structs) {
      my $orf = shift @orf_structs;
      
      my $start = $orf->{start};
      my $stop = $orf->{stop};
      
      my $length = int((abs($start-$stop)+1)/3); 
      my $orient = $orf->{orient};
      my $protein = $orf->{protein};            
      
      ##################################
      # adjust for boundary conditions, since starts and stops run off the ends of the sequences at partial codons
      #################################
      
      # adjust at 3' end
      if ($stop > length($seq)) {
          $stop -= 3;
      }
      if ($start > length($seq)) {
          $start -= 3;
      }
      
      # adjust at 5' end
      if ($stop < 1) {
          $stop += 3;
      }
      if ($start < 1) {
          $start += 3;
      }

      
      if ($length < $MIN_PROT_LENGTH) { next; }
      if ($force_complete and (substr($orf->{protein},0,1) ne 'M'  or substr($orf->{protein},-1) ne '*' ) ) {next;}
      if ($force_start_codon and substr($orf->{protein},0,1) ne 'M' ) {next;}
      

      print "Candidate (len $length): ".Dumper($orf) if $verbose;
      push (@{$canditates{$acc}}, $orf);
      
    }

    if($keep_all_orf){
      my $cpt=1;
      foreach my $orf (@{$canditates{$acc}}){
        #create a new sequence object
        my $new_seqObj  = Bio::Seq->new( '-format' => 'fasta' , -seq => $orf->{protein});
        $new_seqObj->id($acc.".".$cpt);
        $new_seqObj->description($seqObj->description());
        $fasta_out->write_seq($new_seqObj);
        $cpt++;
      }
    }
    else{ # let's keep only the longest 
      my $orf = @{$canditates{$acc}}[0];
      $seqObj->seq($orf->{protein}); #changing the DNA sequence by the corresponding AA sequence is enough
      $fasta_out->write_seq($seqObj);
    }
}

# END
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

#check if reference exists in hash. Deep infinite : hash{a} or hash{a}{b} or hash{a}{b}{c}, etc.
# usage example: exists_keys($hash_omniscient,('level3','cds',$level2_ID)
sub exists_keys {
    my ($hash, $key, @keys) = @_;

    if (ref $hash eq 'HASH' && exists $hash->{$key}) {
        if (@keys) {
            return exists_keys($hash->{$key}, @keys);
        }
        return 1;
    }
    return '';
}

__END__

=head1 NAME

fasta_get_longestORF.pl -
The script take a nucleotide fasta file as input and will extract the longest ORF(s) and translate it(them) in AA.
By default it extracts only the longest ORF even incomplete (missing start or/and stop codon) >= 100 AA.
This script is an adpatation of the TransDecoder.LongestORF tool, adapted to use bioperl.
/!\ Bolean parameter don't expect any value.

=head1 SYNOPSIS

    ./fasta_get_longestORF.pl -f infile.fasta [ -o outfile ]
    ./fasta_get_longestORF.pl -h

=head1 OPTIONS

=over 8

=item B<-f> or B<--fa> or B<--fasta>

Nucleotide fasta file.

=item B<-s> or B<--size_min>

Minimum length of the ORF to be kept in AA (100 by default)

=item B<--ct> or B<--table> or B<--codon>

This option allows specifying the codon table to use - It expects an integer (1 by default = standard)

=item B<--force_start_codon>

This option force to keep the longest ORF that contains a start codon (M). Bolean

=item B<--force_complete>

This option force to keep the longest ORF that contains a start codon (M) and stop codon (*). Bolean

=item B<--keep_all_orf>

This option force to keep all the ORFs that meet the criteria. Bolean

=item B<-v>

Verbose. Useful for debugging purpose. Bolean

=item B<-o> or B<--out> or B<--output> or B<--outfile>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.

=back

=cut
