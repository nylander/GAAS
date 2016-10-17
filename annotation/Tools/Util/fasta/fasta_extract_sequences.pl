#!/usr/bin/env perl

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

fasta_extract_sequences.pl -
This script extract sequence in fasta format from a fasta file. You can extract one fasta sequence providing a sequence name or the name of a file containing a list of sequence name (one by line)

=head1 SYNOPSIS

    ./fasta_extract_sequences.pl -f=infile.fasta -n sequence1 [ -o outfile ]
    ./fasta_extract_sequences.pl --help

=head1 OPTIONS

=over 8

=item B<-f> or B<--fasta> 

Input fasta file.

=item B<-n>, B<--name>

Could be a sequence name to retrieve in the fasta file, or a file containing a list of sequence name (one by line).

=item B<-o> or B<--output>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.

=back

=cut
