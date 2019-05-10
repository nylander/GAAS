#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Pod::Usage;
use IO::File;
use List::MoreUtils qw(uniq);
use File::Basename;
use Bio::Tools::GFF;
use BILS::Handler::GFF3handler qw(:Ok);
use BILS::Handler::GXFhandler qw(:Ok);

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
my $opt_plot = undef;
my $help= 0;

if ( !GetOptions(
    "help|h" => \$help,
    'o|output=s'      => \$opt_output,
    "gff|f=s" => \$gff))

{
    pod2usage( { -message => "Failed to parse command line",
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($help) {
    pod2usage( { -verbose => 1,
                 -exitval => 0,
                 -message => "$header \n" } );
}

if ( ! (defined($gff)) ){
    pod2usage( {
           -message => "$header\nAt least 1 parameter is mandatory:\nInput reference gff file (--gff) \n\n",
           -verbose => 0,
           -exitval => 2 } );
}

#### OUT
my $gffout;
if ($opt_output) {
  my($file, $dirs, $suffix) = fileparse($opt_output, (".gff",".gff1",".gff2",".gff3",".gtf",".gtf1",".gtf2",".gtf3",".txt")); #remove extension 
  if(! $file){print "No output file name provided, just a path...\n";exit;}
  my $path = $dirs.$file;
  open(my $fh, '>', $path.".gff3") or die "Could not open file '$opt_output' $!";
  $gffout= Bio::Tools::GFF->new(-fh => $fh, -gff_version => 3 );
  }
else{
  $gffout = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 3);
}

                #####################
                #     MAIN          #
                #####################



######################
### Parse GFF input #
print "Reading file $gff\n";
my ($hash_omniscient, $hash_mRNAGeneLink) = slurp_gff3_file_JD({ input => $gff
                                                              });
print "Parsing Finished\n";
### END Parse GFF input #
#########################

#check number of level1
my $nbLevel1 = 0;
foreach my $tag_l1 (keys %{$hash_omniscient->{'level1'}}){
  $nbLevel1 += keys %{$hash_omniscient->{'level1'}{$tag_l1}};
}
#chech number of level2
my $nbLevel2 = keys %$hash_mRNAGeneLink;

#Check if we have isoforms
if($nbLevel1 != $nbLevel2){
  

  #create list of level2 where we kept only level2 that have cds and only the longest isoform !
  my $list_id_l2 = get_longest_cds_level2($hash_omniscient);

  # create a new omniscient with only one mRNA isoform per gene
  my $omniscientNew = create_omniscient_from_idlevel2list($hash_omniscient, $hash_mRNAGeneLink, $list_id_l2);

  # print omniscientNew containing only the longest isoform per gene
  print_omniscient($omniscientNew, $gffout);
  print $nbLevel2 - $nbLevel1." isoforms removed ! \n";
}
else{
  print "Nothing to do... this file doesn't contain any isoform !\n";
}

# END STATISTICS #
##################
print "Done\n";




__END__

=head1 NAME

gff3_sp_keep_longest_isoform.pl -
The script take a gff3 file as input. -
The script give basic statistics of a gff file.
Remark: identical feature from level1 or level2 with identical ID will be merged as well as their subsequent features (Level2 or level3).

=head1 SYNOPSIS

    ./gff3_sp_keep_longest_isoform.pl -gff file.gff  [ -o outfile ]
    ./gff3_sp_keep_longest_isoform.pl --help

=head1 OPTIONS

=over 8

=item B<--gff> or B<-f>

Input GFF3 file that will be read (and sorted)

=item B<--output> or B<-o>

File where will be written the result. If no output file is specified, the output will be written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.

=back

=cut
