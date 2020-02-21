#!/usr/bin/env perl

use strict;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use Bio::SeqIO ;
use Bio::DB::Fasta;
use Bio::Tools::GFF;

my $start_run = time();

my $opt_agpfile;
my $opt_fastafile;
my $opt_output;
my $opt_help;

my $header = qq{
########################################################
# NBIS 2018 - Sweden                                   #
# jacques.dainat\@nbis.se                              #
# Please cite NBIS (www.nbis.se) when using this tool. #
########################################################
};

# OPTION MANAGMENT
my @copyARGV=@ARGV;
if ( !GetOptions( 'a|agp=s' => \$opt_agpfile,
                  'f|fa|fasta=s' => \$opt_fastafile,
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

if ( (! (defined($opt_agpfile)) ) or (! (defined($opt_fastafile)) ) ){
    pod2usage( {
           -message => "\nAt least 2 parametes are mandatory:\nInput agp file (-g);  Input fasta file (-f)\n\n".
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

##################################################################
###########################    MAIN    ##########################
##################################################################

######################
### Parse AGP input #
my %hash_agp;
print "Reading file $opt_agpfile\n";
if (open(my $fh, '<:encoding(UTF-8)', $opt_agpfile)) {
  while (my $row = <$fh>) {
    chomp $row;
    my ($object, $object_beg, $object_end, $part_number, $component_type, $component_id_or_gap_length,
        $component_beg_or_gap_type, $component_end_or_linkage, $orientation_or_linkage_evidence ) = split(/\t/, $row);
    push (@{$hash_agp{$object}{$part_number}}, ($object, $object_beg, $object_end, $part_number, $component_type, $component_id_or_gap_length,
        $component_beg_or_gap_type, $component_end_or_linkage, $orientation_or_linkage_evidence));
  }
} else {
  warn "Could not open file '$opt_agpfile' $!";
}

print "Parsing Finished\n";
### END Parse AGP input #
#########################


######################
### READ FASTA input #
my $nbFastaSeq=0;
my $db = Bio::DB::Fasta->new($opt_fastafile);
my @ids      = $db->get_all_primary_ids;
my @newids; foreach my $id (@ids) { $id =~ s/[^[:print:]]//g; push @newids, $id; } # FIX FOR THE CRAZY BUG ADDING NULLBILLION TIMES AT THE END OF MY IDS
my %allIDs; # save ID in lower case to avoid cast problems
foreach my $id (@newids ){$allIDs{lc($id)}=$id;}
### END fASTAS input #
######################

foreach my $object (keys %hash_agp){
  #print $object."\n";
  my $sequence="";
  my $faID="";
  foreach my $part_number (sort {$a <=> $b} keys %{$hash_agp{$object}} ){
    #print $part_number."\n";
    my @feature = @{$hash_agp{$object}{$part_number}};
    $faID = $feature[0];

    if( lc($feature[4]) eq "n" or lc($feature[4]) eq "u" ){ # GAP
      my $gap_length = $feature[5];
      #print "gap length=  $gap_length\n";
      $sequence.= 'N' x $gap_length;
    }
    else{;
      my $peace_sequence = get_sequence($db, \%allIDs, $feature[5], $feature[6],  $feature[7]);
      if($feature[8] eq "+"){
        $sequence.= $peace_sequence;
      }
      elsif($feature[8] eq "-"){
        my $rev_sequence = reverse $peace_sequence;
        $rev_sequence =~ tr/ATCGYRKMDHVBatcgyrkmdhvb/TAGCRYMKHDBVtagcrymkhdbv/;
        $sequence.= $rev_sequence;
      }
      else{
        print "Problem with the strand !!\n";
      }
    }
  }
  #create sequence object
  my $seq  = Bio::Seq->new( '-format' => 'fasta' , -id => $faID , -seq => $sequence);
  $ostream->write_seq($seq);
}

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

sub  get_sequence{
  my  ($db, $allIDs, $seq_id, $start, $end) = @_;

  my $sequence="";
  my $seq_id_correct = undef;
  if( exists_keys($allIDs,(lc($seq_id)) ) ){

    $seq_id_correct = $allIDs{lc($seq_id)};

    $sequence = $db->subseq($seq_id_correct, $start, $end);

    if($sequence eq ""){
      warn "Problem ! no sequence extracted for - $seq_id !\n";  exit;
    }
    if(length($sequence) != ($end-$start+1)){
      my $wholeSeq = $db->subseq($seq_id_correct);
      $wholeSeq = length($wholeSeq);
      warn "Problem ! The size of the sequence extracted ".length($sequence)." is different than the specified span: ".($end-$start+1).".\nThat often occurs when the fasta file does not correspond to the annotation file. Or the index file comes from another fasta file which had the same name and haven't been removed.\n".
           "As last possibility your gff contains location errors (Already encountered for a Maker annotation)\nSupplement information: seq_id=$seq_id ; seq_id_correct=$seq_id_correct ; start=$start ; end=$end ; $seq_id sequence length: $wholeSeq )\n";
    }
  }
  else{
    warn "Problem ! ID $seq_id not found !\n";
  }

  return $sequence;
}

__END__

=head1 NAME

AGP2chromosome.pl -
The script aims to combine contigs from the fasta file in chromosome as described into the AGP file.
AGP version 2 is expected. See https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/ for specification of this format. If you are unsure about the AGP file you are using,
you could check its sanity using the agp validator provided by the NCBI at this address: https://www.ncbi.nlm.nih.gov/projects/genome/assembly/agp/agp_validate.cgi
The result is written to the specified output file, or to STDOUT.


=head1 SYNOPSIS

    ./AGP2chromosome.pl -g=infile.gff -f=infile.fasta  [ -o outfile ]
    ./AGP2chromosome.pl --help

=head1 OPTIONS

=over 8

=item B<--agp> or B<-a>

Input AGP file

=item B<--fasta>, B<--fa>  or B<-f>

Input fasta file.

=item B<-o> or B<--output>

Output GFF file. If no output file is specified, the output will be
written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.

=back

=cut
