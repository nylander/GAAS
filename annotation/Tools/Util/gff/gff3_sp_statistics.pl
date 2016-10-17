#!/usr/bin/env perl

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

#check number of level1
my $nbLevel1 = 0;
foreach my $tag_l1 (keys %{$hash_omniscient->{'level1'}}){
  $nbLevel1 += keys %{$hash_omniscient->{'level1'}{$tag_l1}};
}
#chech number of level2
my $nbLevel2 = keys %$hash_mRNAGeneLink;

##############
# STATISTICS #
my $stat;
if($opt_genomeSize){
  $stat = gff3_statistics($hash_omniscient, $opt_genomeSize);
}else{$stat = gff3_statistics($hash_omniscient);}

#print statistics
foreach my $infoList (@$stat){
  foreach my $info (@$infoList){
    print $out "$info";
  }
  print $out "\n";
}

#Check if we have isoforms
if($nbLevel1 != $nbLevel2){
  print $out "\nApparently we have isoforms : Number level1 features: $nbLevel1 / Number of level2 features $nbLevel2\n";
  print $out "We will proceed to the statistics analysis using only the mRNA with the longest cds\n";

  #create list of level2 where we kept only level2 that have cds and only the longest isoform !
  my $list_id_l2 = get_longest_cds_level2($hash_omniscient);

  # create a new omniscient with only one mRNA isoform per gene
  my $omniscientNew = create_omniscient_from_idlevel2list($hash_omniscient, $hash_mRNAGeneLink, $list_id_l2);
  
  # print stats
  my $stat;
  if($opt_genomeSize){
    $stat = gff3_statistics($omniscientNew, $opt_genomeSize);
  }else{$stat = gff3_statistics($omniscientNew);}

  #print statistics
  foreach my $infoList (@$stat){
    foreach my $info (@$infoList){
      print $out "$info";
    }
    print $out "\n";
  }

}

# END STATISTICS #
##################
print "Bye Bye.\n";
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
sub get_longest_cds_level2{
  my ($hash_omniscient)= @_;

  my @list_id_l2;

  #################
  # == LEVEL 1 == #
  #################
  foreach my $primary_tag_l1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
    foreach my $id_tag_l1 (keys %{$hash_omniscient->{'level1'}{$primary_tag_l1}}){

      #################
      # == LEVEL 2 == #
      #################
      foreach my $primary_tag_l2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
        if ( exists ($hash_omniscient->{'level2'}{$primary_tag_l2}{$id_tag_l1} ) ){

          #check if there is isoforms
          ###########################

          #take only le longest
          if ($#{$hash_omniscient->{'level2'}{$primary_tag_l2}{$id_tag_l1}} > 0){
            my $longestL2 ="";
            my $longestCDSsize = 0;
            foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_l2}{$id_tag_l1}}) {

              my $level2_ID =   lc($feature_level2->_tag_value('ID') ) ;
              if ( exists_keys( $hash_omniscient, ('level3','cds',$level2_ID ) ) ) {

                my $cdsSize=0;
                foreach my $cds ( @{$hash_omniscient->{'level3'}{'cds'}{$level2_ID}} ) { # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
                  $cdsSize += ( $cds->end - $cds->start + 1 );
                }
                if($cdsSize > $longestCDSsize ){
                  $longestL2 = $level2_ID;
                }
              }
            }
            push @list_id_l2,$longestL2; # push id of the longest
          }
          else{ #take it only of cds exits
            my $level2_ID =  lc(@{$hash_omniscient->{'level2'}{$primary_tag_l2}{$id_tag_l1}}[0]->_tag_value('ID')) ;
            if (exists_keys( $hash_omniscient, ('level3','cds', $level2_ID ) ) ){
              push @list_id_l2, $level2_ID; # push the only one existing
            } 
          }
        }
      }
    }
  }

  return \@list_id_l2;
}

__END__

=head1 NAME

gff3_statistics.pl -
The script take a gff3 file as input. -
The script give basic statistics of a gff file.
Remark: identical feature from level1 or level2 with identical ID will be merged as well as their subsequent features (Level2 or level3).

=head1 SYNOPSIS

    ./gff3_statistics.pl -gff file.gff  [ -o outfile ]
    ./gff3_statistics.pl --help

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
