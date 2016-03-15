#!/usr/bin/perl


use Carp;
use Clone 'clone';
use strict;
use Getopt::Long;
use Pod::Usage;
use IO::File;
use Data::Dumper;
use List::MoreUtils qw(uniq);
use Bio::Tools::GFF;
use BILS::Handler::GXFhandler qw(:Ok);
use BILS::Handler::GFF3handler qw(:Ok);
use BILS::GFF3::Statistics qw(:Ok);

my $header = qq{
########################################################
# BILS 2015 - Sweden                                   #  
# jacques.dainat\@bils.se                               #
# Please cite BILS (www.bils.se) when using this tool. #
########################################################
};

my $gff = undef;
my $opt_output = undef;
my $opt_genomeSize = undef;
my $help= 0;

if ( !GetOptions(
    "help|h" => \$help,
    'o|output=s'      => \$opt_output,
    'g|gs=s' => \$opt_genomeSize,
    "gff|f=s" => \$gff))

{
    pod2usage( { -message => '$header Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($help) {
    pod2usage( { -verbose => 2,
                 -exitval => 2,
                 -message => "$header \n" } );
}
 
if ( ! (defined($gff)) ){
    pod2usage( {
           -message => "$header\nAt least 1 parameter is mandatory:\nInput reference gff file (--gff) \n\n",
           -verbose => 0,
           -exitval => 1 } );
}


#### IN / OUT
my $out = IO::File->new();
if ($opt_output) {

  if (-f $opt_output){
      print "Cannot create a directory with the name $opt_output because a file with this name already exists.\n";exit();
  }
  if (-d $opt_output){
      print "The output directory choosen already exists. Please geve me another Name.\n";exit();
  }

  open($out, '>', $opt_output) or die "Could not open file '$opt_output' $!";
  }
else{
  $out->fdopen( fileno(STDOUT), 'w' );
}
                #####################
                #     MAIN          #
                #####################



######################
### Parse GFF input #
print "Reading file $gff\n";
my ($hash_omniscient, $hash_mRNAGeneLink) = BILS::Handler::GFF3handler->slurp_gff3_file_JD($gff);
print "Parsing Finished\n";
### END Parse GFF input #
#########################


##############
# STATISTICS #
my $stat;
if($opt_genomeSize){
  $stat = gff3_statistics($hash_omniscient, $opt_genomeSize);
}else{$stat = gff3_statistics($hash_omniscient);}

#print statistics
foreach my $info (@$stat){
  print $out "$info";
}
# END STATISTICS #
##################

__END__

=head1 NAME

gff3_checkOmniscient.pl -
The script take a gff3 file as input. -
The script give basic statistics of a gff file. 
Remark: identical feature from level1 or level2 with identical ID will be merged as well as their subsequent features (Level2 or level3).

=head1 SYNOPSIS

    ./gff3_checkOmniscient.pl -gff file.gff  [ -o outfile ]
    ./gff3_checkOmniscient.pl --help

=head1 OPTIONS

=over 8

=item B<--gff> or B<-f>

Input GFF3 file that will be read (and sorted)

=item B<--gs> or B<-g>

This option inform about the genome size in oder to compute more statistics. You can give the size in Nucleotide or directly the fasta file.

=item B<--output> or B<-o>

File where will be written the result. If no output file is specified, the output will be written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.

=back

=cut