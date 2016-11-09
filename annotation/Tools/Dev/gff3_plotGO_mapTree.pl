#!/usr/bin/env perl

#####################################################################
# maker_checkFusionSplitBetweenTwoBuilds v1 - Jacques Dainat 10/2014 #
#####################################################################

use strict;
use warnings;
use POSIX qw(strftime);
use Data::Dumper;
use Carp;
use Getopt::Long;
use IO::File;
use Pod::Usage;
use Statistics::R;
use BILS::Handler::GXFhandler qw(:Ok);
use BILS::Handler::GFF3handler qw(:Ok);
use Bio::OntologyIO;
use Bio::Tools::GFF;

my $header = qq{
########################################################
# BILS 2015 - Sweden                                   #  
# jacques.dainat\@bils.se                               #
# Please cite BILS (www.bils.se) when using this tool. #
########################################################
};


my $opt_reffile;
my $opt_obo=undef;
my $opt_output=undef;
my $opt_help = 0;
my $DefaultUTRnb=5;

my @copyARGV=@ARGV;
if ( !GetOptions( 'f|gff|ref|reffile=s' => \$opt_reffile,
                  'obo=s' => \$opt_obo,
                  'o|out|output=s' => \$opt_output,            
                  'h|help!'         => \$opt_help ) )
{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

if ($opt_help) {
    pod2usage( { -verbose => 2,
                 -exitval => 2,
                 -message => "$header\n" } );
}

if ( ! defined($opt_reffile ) or ! $opt_output or ! $opt_obo) {
    pod2usage( {
           -message => "$header\nMust specify at least 2 parameters:\nReference data gff3 file (--gff)\nOne UTR option (3, 5 , both, plot)",
           -verbose => 0,
           -exitval => 1 } );
}

# #####################
# # Manage Input File #
# #####################
 my $ref_istream = Bio::Tools::GFF->new(-file => $opt_reffile, -gff_version => 3 ) or
      croak(sprintf( "Can not open '%s' for writing %s", $opt_reffile, $! ));

# #########################
# # END Manage Input File #
# #########################
# #######################
# # START Manage Option #
# #######################
if (-f $opt_output){
    print "The output directory choosen already exists. Please geve me another Name.\n";exit();
}

my  $ostreamReport = \*STDOUT or die ( sprintf( "Can not open '%s' for writing %s", "STDOUT", $! ));

my $string1 = strftime "%m/%d/%Y at %Hh%Mm%Ss", localtime;
$string1 .= "\n\nusage: $0 @copyARGV\n\n";

print $ostreamReport $string1;

# #####################################
# # END Manage OPTION  
# #####################################
# print "parse obo file\n";
# my $parser = Bio::OntologyIO->new
#     ( -format       => "obo",
#       -file        =>  $opt_obo);

# my %hash_onto;
# while(my $ont = $parser->next_ontology()) {
# 	#print Dumper($ont);
# 	my $ontoName= $ont->name();
# 	print "read ontology ",$ont->name()," with ",
# 	scalar($ont->get_root_terms)," root terms, and ",
#     scalar($ont->get_all_terms)," total terms, and ",
#     scalar($ont->get_leaf_terms)," leaf terms\n";
#     $hash_onto{$ontoName} = $ont;
# }
# my $hash_bioprocess = $hash_onto{'biological_process'};
# my @term = $hash_bioprocess->find_terms(-identifier => "GO:0071428");
# #print Dumper($term);
#   print $term[0]->identifier(), "\n";
#   print $term[0]->name(), "\n";
# #  print $term[0]->description(), "\n";
# #  print $term->definition(), "\n";
# #  print $term->is_obsolete(), "\n";
# #  print $term->comment(), "\n";
# exit;
                                                      #######################
                                                      #        MAIN         #
#                     >>>>>>>>>>>>>>>>>>>>>>>>>       #######################       <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


# #################################
# # Manage Ouput Directory / File #
# #################################


######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) = BILS::Handler::GXFhandler->slurp_gff3_file_JD($opt_reffile);
print("Parsing Finished\n\n");
### END Parse GFF input #
#########################

###########################
# get GO terms information
###########################
my %GOdist;
foreach my $tag_l1 (keys %{$hash_omniscient->{'level1'}}){ 
  foreach my $id_l1 (keys %{$hash_omniscient->{'level1'}{$tag_l1}}){
      
    #################
    # == LEVEL 2 == #
    #################
    foreach my $tag_l2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...        
      if ( exists ($hash_omniscient->{'level2'}{$tag_l2}{$id_l1} ) ){
        foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}}) {
          if($feature_level2->has_tag('Ontology_term')){
            my @GOterms=$feature_level2->get_tag_values('Ontology_term');
            foreach my $GOterm (@GOterms){
            	$GOdist{$GOterm}++;
            }
          }
          else{
            $GOdist{"none"}++;
          }
        }
      }
    }
  }
}


#####################
# Plot distribution 




my $txtFile;
my $outPlot;

$txtFile = $opt_output.".txt";
$outPlot = $opt_output.".pdf";

#print file thtat will be read by R
open(FH, ">".$txtFile) || die "Erreur E/S:$!\n";
my $firstLine="yes";
foreach my $GO_type (keys %GOdist) {
  if($firstLine){
     print FH $GO_type."\t".0;
     $firstLine=undef;
  }else{
     print FH "\n".$GO_type."\t".0;
  }
}
close FH;
exit;
my $R = Statistics::R->new() or die "Problem with R : $!\n";

#R command
#$R->run(q`install.packages( "treemap" );`);
#$R->run(q`library(treemap) `);

 $R->send(
     qq`
    #install library if needed
    if(!require(treemap)){
     	chooseCRANmirror()
     	options("repos")
  		install.packages("treemap")
  		library(treemap)
		}
	library(treemap)`
         );



#     listValues=as.matrix(read.table("$txtFile", sep="\t", he=F)) ##///!!!\\\\\
#     legendToDisplay=paste("Number of value used : ",length(listValues))
#     listValueMoreThan <- listValues[listValues[,1]>5,]

#     pdf("$outPlot")
#     plot(listValues[,2]~listValues[,1], xlab="Contig size", ylab="Frequency", main="Size distribution of $utr_type")
#     dev.off()
    
#     pdf("$outPlotOver")
#     plot(listValueMoreThan[,2]~listValueMoreThan[,1], xlab="Contig size", ylab="Frequency", main="Size distribution of $utr_type over 5")
#     dev.off()`
#         );

# Close the bridge
$R->stopR();

# Delete temporary file
unlink "$txtFile";





#######################################################################################################################
        ####################
         #     methods    #
          ################
           ##############
            ############
             ##########
              ########
               ######
                ####
                 ##

__DATA__ 


__END__


=head1 NAME

maker_manageUTR.pl - Detect the genes containing too much UTR's exon according to a choosen threshold.
If no UTR option (3, 5, 3 and 5, both) is given the threshold will be not used. 
option 3 and 5 together is different of "both". In the first case the gene is discarded if either the 3' or the 5' UTR contains more exon than the threshold given.
In the second case, will be discarded only the genes where the addition of UTR's exon of both side is over the threshold given.

=head1 SYNOPSIS

    ./maker_manageUTR.pl --ref=infile --three --five -p --out=outFile 
    ./maker_manageUTR.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>, B<--ref>, B<--reffile> or B<-f>

Input GFF3 file correponding to gene build.

=item B<-n>, B<--nb> or B<--number>

Threshold of exon's number of the UTR. Over or equal to this threshold, the UTR will be discarded. Default value is 5.

=item B<-3>, B<--three> or B<--tree_prime_utr>

The threshold of the option <n> will be applied on the 3'UTR.

=item B<-5>, B<--five> or B<--five_prime_utr>

The threshold of the option <n> will be applied on the 5'UTR.

=item B<-b>, B<--both> or B<--bs>

The threshold of the option <n> will be applied on genes where the number of UTR exon (3' and 5' additioned) is over it.

=item  B<--p>, B<--plot> or B<-o>

Allows to create an histogram in pdf of UTR sizes distribution.

=item  B<--out>, B<--output> or B<-o>

Output gff3 file where the gene incriminated will be write.

=item B<--help> or B<-h>

Display this helpful text.

=back

=cut
