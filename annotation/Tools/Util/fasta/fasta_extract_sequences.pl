#!/usr/bin/perl

use strict;
use Pod::Usage;
use Getopt::Long;
use Bio::SeqIO ;
use Bio::DB::Fasta;
use Data::Dumper;

my $start_run = time();

my $opt_fastafile;
my $opt_output;
my $opt_help = 0;
my $opt_name = undef; 

my $header = qq{
########################################################
# BILS 2016 - Sweden                                   #  
# jacques.dainat\@bils.se                               #
# Please cite BILS (www.bils.se) when using this tool. #
########################################################
};

# OPTION MANAGMENT
my @copyARGV=@ARGV;
if ( !GetOptions( 'f|fa|fasta=s' => \$opt_fastafile,
                  'n|name=s' => \$opt_name,
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
 
if ( (! (defined($opt_name)) ) or (! (defined($opt_fastafile)) ) ){
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
#### read fasta file and save info in memory
######################
my $db = Bio::DB::Fasta->new($opt_fastafile);
print ("Genome fasta parsed\n");

my @list_seq_result=();
if(-f $opt_name){
  print "It's a file you gave me... I will read line by line to take in account all the names you proveded\n";
  print "To be implemented .... sorry\n";
}
else{
  if($opt_name =~ m/^>/){
    print "I remove the chevron (\">\") !\n";
    $opt_name=substr($opt_name, 1 , length($opt_name));
    print $opt_name."\n";
  }
  if($db->seq($opt_name)){
    #create sequence object
    my $seq_obj  = Bio::Seq->new( '-format' => 'fasta' , -seq => $db->seq($opt_name), -display_id => $opt_name);
    push @list_seq_result, $seq_obj;
  }
}

if (! @list_seq_result){
  print "Nothing found !\n";
}
else{
  foreach my $seq_obj (@list_seq_result){
    $ostream->write_seq($seq_obj);  
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

gff3_extract_cds_sequences.pl -
This script extract sequence in feasta format from gff file. You can extract the fasta of any kind of feature define by the 2th column in the gff file.
The result is written to the specified output file, or to STDOUT.

The Header are formated like that:
>mRNA_ID|gene_ID|Name=NAME|Seq_id=Chromosome_ID|type=cds|5'extra=VALUE

/!\ mRNA_ID not displayed when extracting gene. 
Name is optional and will be written only if the Name attribute exists in th gff.
type will be the feature type extracted. 
5'extra or 3'extra is otpional, according to the use of the upstream and downstream option.

=head1 SYNOPSIS

    ./gff3_extract_cds_sequences.pl -g=infile.gff -f=infile.fasta  [ -o outfile ]
    ./gff3_extract_cds_sequences.pl --help

=head1 OPTIONS

=over 8

=item B<-g>, B<--gff> or B<-ref>

Input GFF3 file that will be read (and sorted)

=item B<-f> or B<--fasta> 

Input fasta file.

=item B<-t> 

Define the feature you want to extract the sequnce from. By deafault it's 'cds'. Most common choice are: gene,mrna,exon,cds,trna,three_prime_utr,five_prime_utr.

=item B<-p>, B<--protein> or B<--aa>

Will translate the extracted sequence in Amino acid.

=item B<-e> or B<--ext>

This option called "extremities" allows dealing with multifeature like cds or exon, to extract the full sequence from start extremity to the end extremity, i.e with introns.
Use of that option with exon will give the same result as extract the mrna sequence.
Use of that option on cds will give the cdna wihtout the untraslated sequences.

=item B<-u>, B<--up>, B<-5>, B<--five> or B<-upstream>

Integer. It will take that number of nucleotide in more at the 5' extremity. Option "e" must be activated to use this option (Why ? to avoid to extract intronic/overlaping sequence in case of feature spread over several locations (exon,cds,utrs)).


=item B<-d>, B<--do>, B<-3>, B<--three>, B<-down> or B<-downstream>

Integer. It will take that number of nucleotide in more at the 3' extremity. Option "e" must be activated to use this option (Why ? to avoid to extract intronic/overlaping sequence in case of feature spread over several locations (exon,cds,utrs)).

=item B<--cdna>

This extract the cdna sequence (i.e transcribed sequence (devoid of introns, but containing untranslated exons)). It correspond to extract the exons sequences.

=item B<-o> or B<--output>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.

=back

=cut
