#!/usr/bin/env perl

###
# Implement case insensitive
###
use strict;
use Pod::Usage;
use Getopt::Long;
use Bio::SeqIO ;
use Bio::DB::Fasta;
use GAAS::GAAS;

my $header = get_gaas_header();
my $start_run = time();

my $col = undef;
my $lineToAvoid=undef;
my $separator=undef;
my $opt_fastafile;
my $opt_output;
my $opt_help = 0;
my $opt_name = undef;

# OPTION MANAGMENT
my @copyARGV=@ARGV;
if ( !GetOptions( 'f|fa|fasta=s' => \$opt_fastafile,
                  "line=i" => \$lineToAvoid,
                  "col=i" => \$col,
                  "s=s" =>\$separator,
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
    pod2usage( { -verbose => 99,
                 -exitval => 0,
                 -message => "$header\n" } );
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


###########################
### Extract sequence ID(s)
my %list_of_ID;
my $nbID=0;
# Case it's a file
if (-f $opt_name){
  #Manage column with the ID
  if (! defined $col){
    $col=0;
  }
  else{$col=$col -1 ;}

  #Manage line to avoid
  if (! defined $lineToAvoid){
    $lineToAvoid=0;
  }

  print "It's a file you gave me... I will read line by line to take in account the ID from the column ".($col+1)." of the file $opt_name\n";
  print "The first $lineToAvoid lines will be ignored.\n";
  my $ID_list  = IO::File->new("<".$opt_name);


  my $cpt_line=0;
  while ( <$ID_list> ) {
    $cpt_line++;

    if($cpt_line > $lineToAvoid){

      chomp;
      if(! $_ =~ /^\s*$/){

      my @cols;
        if (! $separator){
          @cols = split /\s/, $_;
      }
      else{
        @cols = split /$separator/, $_;
      }
        my $id = $cols[$col];
        $id =~ s/[^[:print:]]+//g;
        #print $id."\n";
        if($id =~ m/^>/){
          print "I remove the chevron (\">\") !\n";
          $id=substr($id, 1 , length($id));
        }
        $list_of_ID{$id}++;
        $nbID++;
      }
    }
  }
  print "$nbID ID found.\n";
}
else{
  if($opt_name =~ m/^>/){
    print "I remove the chevron (\">\") !\n";
    $opt_name=substr($opt_name, 1 , length($opt_name));
    print $opt_name."\n";
  }
  $list_of_ID{$opt_name}++;
}

##########################
#Now extract the sequences
my @list_seq_result=();
foreach my $ID (keys %list_of_ID){

  if($db->seq($ID)){
    #create sequence object
    my $seq_obj = $db->get_Seq_by_id($ID);

    push @list_seq_result, $seq_obj;
  }
  else{
    print "<$ID> not found into the $opt_fastafile fasta file !\n";
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

gaas_fasta_extract_sequence_from_id.pl

=head1 DESCRIPTION

This script extract sequence in fasta format from a fasta file. You can extract one fasta sequence providing a sequence name or the name of a file containing a list of sequence name (one by line)

=head1 SYNOPSIS

    gaas_fasta_extract_sequence_from_id.pl -f=infile.fasta -n sequenceID [ -o outfile ]
    gaas_fasta_extract_sequence_from_id.pl --help

=head1 OPTIONS

=over 8

=item B<-f> or B<--fasta>

Input fasta file.

=item B<-n>, B<--name>

Could be a sequence name to retrieve in the fasta file, or a file containing a list of sequence name (one by line).

=item B<--line>

Integer, number of line to avoid. Allow to avoid headers.

=item B<--col>

column containing the ID. By default the first column is considered.

=item B<-q>

Field separator, by default un-printable character are use as separator (\s). You can define the one you wnat with this option.

=item B<-o> or B<--output>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

=item B<-h> or B<--help>

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
