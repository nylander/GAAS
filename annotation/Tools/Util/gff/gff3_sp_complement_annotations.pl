#!/usr/bin/env perl


use Carp;
use Clone 'clone';
use strict;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use List::MoreUtils qw(uniq);
use Bio::Tools::GFF;
use BILS::Handler::GFF3handler qw(:Ok);
use BILS::Handler::GXFhandler qw(:Ok);

my $header = qq{
########################################################
# BILS 2016 - Sweden                                   #  
# jacques.dainat\@bils.se                               #
# Please cite BILS (www.bils.se) when using this tool. #
########################################################
};

my $outfile = undef;
my @opt_files;
my $file2 = undef;
my $help= 0;

if ( !GetOptions(
    "help|h" => \$help,
    "gff|f=s" => \@opt_files,
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

if ( ! @opt_files or (@opt_files and ($#opt_files < 1) ) ){
    pod2usage( {
           -message => "\nAt least 2 files are mandatory:\n --gff file1 --gff file2\n\n",
           -verbose => 0,
           -exitval => 2 } );
}

######################
# Manage output file #
my $gffout;
if ($outfile) {
  $outfile=~ s/.gff//g;
open(my $fh, '>', $outfile.".gff") or die "Could not open file '$outfile' $!";
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

my $file1 = shift @opt_files;
my ($hash_omniscient, $hash_mRNAGeneLink) = BILS::Handler::GXFhandler->slurp_gff3_file_JD($file1);
print ("$file1 GFF3 file parsed\n");
info_omniscient($hash_omniscient);

#Add the features of the other file in the first omniscient. It takes care of name to not have duplicates
foreach my $next_file (@opt_files){
  my ($hash_omniscient2, $hash_mRNAGeneLink2) = BILS::Handler::GXFhandler->slurp_gff3_file_JD($next_file);
  print ("$next_file GFF3 file parsed\n");
  info_omniscient($hash_omniscient2);

  ################################
  # First rename ID to be sure to not add feature with ID already used
  rename_ID_existing_in_omniscient($hash_omniscient, $hash_omniscient2);
  print ("\n$next_file IDs checked and fixed.\n");


  # Quick stat hash before complement
  my %quick_stat1;
  foreach  my $tag (keys %{$hash_omniscient->{'level1'}}) {
    my $nb_tag = keys %{$hash_omniscient->{'level1'}{$tag}};
    $quick_stat1{$tag} = $nb_tag;
  }

  ####### COMPLEMENT #######
  complement_omniscients($hash_omniscient, $hash_omniscient2);
  print ("\nComplement done !\n");


 #RESUME COMPLEMENT
  my $complemented=undef;
  # Quick stat hash after complement
  my %quick_stat2;
  foreach  my $tag (keys %{$hash_omniscient->{'level1'}}) {
    my $nb_tag = keys %{$hash_omniscient->{'level1'}{$tag}};
    $quick_stat2{$tag} = $nb_tag;
  }

  #About tag from hash1 added which exist in hash2
  foreach my $tag (keys %quick_stat1){
    if ($quick_stat1{$tag} != $quick_stat2{$tag} ){
      print "We added ".($quick_stat2{$tag}-$quick_stat1{$tag})." $tag(s)\n";
      $complemented=1;
    }
  }
  #About tag from hash2 added which dont exist in hash1
  foreach my $tag (keys %quick_stat2){
    if (! exists $quick_stat1{$tag} ){
      print "We added".$quick_stat2{$tag}."\n";
      $complemented=1;
    }
  }
  #If nothing added
  if(! $complemented){
    print "\nNothing has been added\n";
  }
  else{
    print "\nNow the data contains:\n";
    info_omniscient($hash_omniscient);
  }
}

########
# Print results
print_omniscient($hash_omniscient, $gffout);  

__END__

=head1 NAME
 
gff3_sp_complement_annotations.pl - 
This script merge different gff annotation files in gff format in one. It takes care of duplicated names, and merge overlaping genes together.

=head1 SYNOPSIS

    ./gff3_sp_complement_annotations.pl --gff=infile1 --gff=infile2 --out=outFile 
    ./gff3_sp_complement_annotations.pl --help

=head1 OPTIONS

=over 8

=item B<--gff> or B<-f>

Input GFF3 file(s). You can specify as much file you want like so: -f file1 -f file2 -f file3

=item  B<--out>, B<--output> or B<-o>

Output gff3 file where the gene incriminated will be write.

=item B<--help> or B<-h>

Display this helpful text.

=back

=cut
