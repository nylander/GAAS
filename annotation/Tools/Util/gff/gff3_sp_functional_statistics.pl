#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Pod::Usage;
use Statistics::R;
use IO::File;
use Data::Dumper;
use List::MoreUtils qw(uniq);
use Bio::Tools::GFF;
use BILS::Handler::GFF3handler qw(:Ok);
use BILS::Handler::GXFhandler qw(:Ok);
use BILS::GFF3::Statistics qw(:Ok);

my $header = qq{
########################################################
# BILS 2019 - Sweden                                   #
# jacques.dainat\@nbis.se                               #
# Please cite NBIS (www.nbis.se) when using this tool. #
########################################################
};

my $gff = undef;
my $opt_output = undef;
my $opt_genomeSize = undef;
my $help= 0;

if ( !GetOptions(
    "help|h" => \$help,
    'g|gs=s' => \$opt_genomeSize,
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
my ($hash_omniscient, $hash_mRNAGeneLink) = BILS::Handler::GXFhandler->slurp_gff3_file_JD($gff, undef, undef, undef);
print "Parsing Finished\n";
### END Parse GFF input #
#########################

###############################################################
### Print Statistics structural first
###############################################################
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
my $distri;
if($opt_genomeSize){
  ($stat, $distri) = gff3_statistics($hash_omniscient, $opt_genomeSize);
}
else{
  ($stat, $distri) = gff3_statistics($hash_omniscient);
}

#print statistics
foreach my $infoList (@$stat){
  foreach my $info (@$infoList){
    print $out "$info";
  }
  print $out "\n";
}

# #Check if we have isoforms
# if($nbLevel1 != $nbLevel2){

#   print $out "\nApparently we have isoforms : Number level1 features: $nbLevel1 / Number of level2 features $nbLevel2\n";
#   print $out "We will proceed to the statistics analysis using only the mRNA with the longest cds\n";

#   #create list of level2 where we kept only level2 that have cds and only the longest isoform !
#   my $list_id_l2 = get_longest_cds_level2($hash_omniscient);

#   # create a new omniscient with only one mRNA isoform per gene
#   my $omniscientNew = create_omniscient_from_idlevel2list($hash_omniscient, $hash_mRNAGeneLink, $list_id_l2);
  
#   # print stats
#   my $stat;
#   my $distri;
#   if($opt_genomeSize){
#     ($stat, $distri) = gff3_statistics($omniscientNew, $opt_genomeSize);
#   }else{ 
#     ($stat, $distri) = gff3_statistics($omniscientNew);
#   }

#   #print statistics
#   foreach my $infoList (@$stat){
#     foreach my $info (@$infoList){
#       print $out "$info";
#     }
#     print $out "\n";
#   }
# }

###############################################################
### Print Statistics function 
###############################################################:
my %names_l1=undef;
my $name_l1_nb=undef;
my %names_l2=undef;
my $name_l2_nb=undef;
my %products=undef;
my $product_l2_nb=undef;
my %descriptions=undef;
my $description_l2_nb=undef;
my %ontology_terms=undef;
my $ontology_term_l2_nb=undef;

my %DB_omni_mrna=undef ;
my %DB_omni_gene=undef ;

my $nbmRNAwithFunction=0;
my $nbGeneWithFunction=0;

  #################
  # == LEVEL 1 == #
  #################
foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
  foreach my $gene_id_tag_key (keys %{$hash_omniscient->{'level1'}{$primary_tag_key_level1}}){
    my $l1_have_function=undef;
    my $gene_feature=$hash_omniscient->{'level1'}{$primary_tag_key_level1}{$gene_id_tag_key};
    my $id_gene=$gene_feature->_tag_value('ID');

    #Check For NAME 
    if($gene_feature->has_tag('Name') ){
      my $value = $gene_feature->_tag_value('Name');
      $names_l1{$value}++;
      $name_l2_nb++
      #print "l1 has tag name with value:".$value."\n";
    }



    foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
      if ( exists_keys( $hash_omniscient, ('level2', $primary_tag_key_level2, $gene_id_tag_key) ) ){
        foreach my $level2_feature ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$gene_id_tag_key}}) {
          my $l2_have_function=undef;
          my $id_mrna=$level2_feature->_tag_value('ID');

          #Check For NAME 
          if($level2_feature->has_tag('Name') ){
            my $value = $level2_feature->_tag_value('Name');
            $names_l2{$value}++;
            $name_l2_nb++;
            #print "l2 has tag name with value:".$value."\n";
          }

          #Check For product 
          if($level2_feature->has_tag('product') ){
            my $value = $level2_feature->_tag_value('product');
            $products{$value}++;
            $product_l2_nb++;
            $l2_have_function=1;
            #print "l2 has tag product with value:".$value."\n";
          }
     

          #Check For descritpion 
          if($level2_feature->has_tag('description') ){
            my $value = $level2_feature->_tag_value('description');
            $descriptions{$value}++;
            $description_l2_nb++;
            $l2_have_function=1;
            #print "l2 has tag descritpion with value:".$value."\n";
          }

          #Check For Ontology_term 
          if($level2_feature->has_tag('Ontology_term') ){
            my @value = $level2_feature->_tag_value('Ontology_term');
            $ontology_terms{$value}++;
            $ontology_term_l2_nb++;
            $l2_have_function=1;
            $DB_omni_mrna{'Ontology_term'}{$id_mrna}++;
            $DB_omni_gene{'Ontology_term'}{$id_gene}++;
            #print "l2 has tag ontology_term with value:".$value."\n";
          }

          #Check For Dbxref 
          if($level2_feature->has_tag('Dbxref') ){
            my @values = $level2_feature->get_tag_values('Dbxref');
            foreach my $tuple (@values){
              my ($type,$value) = split /:/,$tuple;
              print $type." ".$value."\n";
              if($type == ""){print "stop";exit;}
              $DB_omni_mrna{$type}{$id_mrna}++;
              $DB_omni_gene{$type}{$id_gene}++;
            }
            $l2_have_function=1;
          }

          if($l2_have_function){
            $nbmRNAwithFunction++;
            $l1_have_function=1;
          }
        }
      }
    }
    if($l1_have_function){
      $nbGeneWithFunction++;
    }
  }
}

my $nbmRNAwithoutFunction=0;
my $nbGeneWithoutFunction=0;

my $listOfFunction;
foreach my $funct (sort keys %DB_omni_mrna){
  $listOfFunction.="$funct,";
}
chop $listOfFunction;

# NOW summerize
my $stringPrint=undef;

my $lineB=       "___________________________________________________________________________________________________";
$stringPrint .= " ".$lineB."\n";
$stringPrint .= "|          | Nb term linked to mRNA | Nb mRNA with term | Nb gene with term |\n";
$stringPrint .= "|".$lineB."|\n";

foreach my $type (keys %DB_omni_mrna){
    my $total_term_mRNA=0;
    foreach my $id_l2 (keys %{$DB_omni_mrna{$type}} ){
      $total_term_mRNA+=$DB_omni_mrna{$type}{$id_l2};
    }
    my $nbmRNA_with_term = keys %{$DB_omni_mrna{$type}};
    my $nbGenewith_term = keys %{$DB_omni_gene{$type}};

    my $mRNA_type =0; #keys %{$mRNAAssociatedToTerm{$type}};
    my $gene_type =0; #keys %{$GeneAssociatedToTerm{$type}};
    $stringPrint .= "|".sizedPrint(" $type",10)."|".sizedPrint($total_term_mRNA,25)."|".sizedPrint($nbmRNA_with_term,25)."|".sizedPrint($nbGenewith_term,25)."|\n|".$lineB."|\n";
  }


$stringPrint .= "nb mRNA without Functional annotation ($listOfFunction) = $nbmRNAwithoutFunction;\n".
                  "nb mRNA with Functional annotation ($listOfFunction) = $nbmRNAwithFunction\n".
                  "nb gene without Functional annotation ($listOfFunction) = $nbGeneWithoutFunction\n".
                  "nb gene with Functional annotation ($listOfFunction) = $nbGeneWithFunction\n"; 


print $out $stringPrint;
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


sub sizedPrint{
  my ($term,$size) = @_;
  my $result; 
  my $sizeTerm = ($term) ? length($term) : 0;
  if ($sizeTerm > $size ){
    $result=substr($term, 0,$size);
    return $result;
  }
  else{
    my $nbBlanc=$size-$sizeTerm;
    $result=$term;
    for (my $i = 0; $i < $nbBlanc; $i++){
      $result.=" ";
    }
    return $result;
  }
}


__END__

########################
# Manage Interpro File #
if (defined $opt_InterproFile){
  parse_interpro_tsv($streamInter,$opt_InterproFile);
  
  # create streamOutput
  if($opt_output){
    foreach my $type (keys %functionData){
      my $ostreamFunct = IO::File->new(); 
      $ostreamFunct->open( $opt_output."/$type.txt", 'w' ) or
          croak(
              sprintf( "Can not open '%s' for writing %s", $opt_output."/$type.txt", $! )
          );
      $functionStreamOutput{$type}=$ostreamFunct;
    }
  }
}
# END MANAGE FUNCTIONAL INPUT FILE #
####################################

##############################
# print FUNCTIONAL INFORMATION

# first table name\tfunction
if($opt_output){
  foreach my $function_type (keys %functionOutput){
    my $streamOutput=$functionStreamOutput{$function_type};
    foreach my $ID (keys %{$functionOutput{$function_type}}){

        print $streamOutput  $ID."\t".$functionOutput{$function_type}{$ID}."\n";

    }
  }
}


# NOW summerize
$stringPrint =""; # reinitialise (use at the beginning)
if ($opt_InterproFile){
  #print INFO
  my $lineB=       "___________________________________________________________________________________________________";
  $stringPrint .= " ".$lineB."\n";
  $stringPrint .= "|          | Nb Total term | Nb mRNA with term  | Nb mRNA updated by term | Nb gene updated by term |\n";
  $stringPrint .= "|          | in Annie File |   in Annie File    | in our annotation file  | in our annotation file  |\n";
  $stringPrint .= "|".$lineB."|\n";

  foreach my $type (keys %functionData){
    my $total_type = $TotalTerm{$type};
    my $mRNA_type_Annie = $functionDataAdded{$type};
    my $mRNA_type = keys %{$mRNAAssociatedToTerm{$type}};
    my $gene_type = keys %{$GeneAssociatedToTerm{$type}};
    $stringPrint .= "|".sizedPrint(" $type",10)."|".sizedPrint($total_type,15)."|".sizedPrint($mRNA_type_Annie,20)."|".sizedPrint($mRNA_type,25)."|".sizedPrint($gene_type,25)."|\n|".$lineB."|\n";
  }

  #RESUME TOTAL OF FUNCTION ATTACHED
  my $listOfFunction;
  foreach my $funct (keys %functionData){
    $listOfFunction.="$funct,";
  }
  chop $listOfFunction;
  my $nbGeneWithoutFunction= keys %geneWithoutFunction;
  my $nbGeneWithFunction= keys %geneWithFunction;
  $stringPrint .= "nb mRNA without Functional annotation ($listOfFunction) = $nbmRNAwithoutFunction\n".
                  "nb mRNA with Functional annotation ($listOfFunction) = $nbmRNAwithFunction\n".
                  "nb gene without Functional annotation ($listOfFunction) = $nbGeneWithoutFunction\n".
                  "nb gene with Functional annotation ($listOfFunction) = $nbGeneWithFunction\n"; 
}


  #Lets keep track the duplicated names
  if($opt_output){
    my $duplicatedNameOut=IO::File->new(">".$opt_output."/duplicatedNameFromBlast.txt" );
    foreach my $name (sort { $duplicateNameGiven{$b} <=> $duplicateNameGiven{$a} } keys %duplicateNameGiven){
      print $duplicatedNameOut "$name\t".($duplicateNameGiven{$name}+1)."\n";
    }
  }

=head1 NAME

gff3_statistics.pl -
The script take a gff3 file as input. -
The script give basic statistics of a gff file.
Remark: identical feature from level1 or level2 with identical ID will be merged as well as their subsequent features (Level2 or level3).

=head1 SYNOPSIS

    ./gff3_statistics.pl --gff file.gff  [ -o outfile ]
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
