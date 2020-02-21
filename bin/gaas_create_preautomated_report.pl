#!/usr/bin/env perl

#LIBRARIES
use strict;
use warnings;
use Data::Dumper;
use YAML::XS 'LoadFile';
use Getopt::Long;
use Pod::Usage;
use IO::File;

#PARAMETERS
my $input_file;
my $output_file;

#VARIABLES
my $config_file;
my $config;
my $fasta_file;
my $repeat_file;
my $gene_file;
my $cegma_file;
my $function_file;
my $species;
my $strain;
my $otherRNA_file;
my $directory ='.';

#a little hash
my %metrics = (
    'species',
    'strain',
    'fasta' => {
      'seq_number' => {},
      'nuc_number' => {},
      'N_number' => {},
      'N_region_number' => {},
      'GC_content' => {},
      'N50' => {}
    },
    'cegma' => {
      'prot_complete_number' => {},
      'percent_completeness' => {},
      'complete_total' => {},
      'complete_average' => {},
      'percent_complete_orthologs' => {},
      'prot_partial_number' => {},
      'percent_partial' => {},
      'partial_total' => {},
      'partial_average' => {},
      'percent_partial_orthologs' => {}
    },
    'repeats' => {
      'number'=> {},
      'total_size'=> {},
      'mean_size'=> {},
      'percent_genome'=> {}
    },
    'gene_build' => {
      'prot_coding_gene_nb' => {},
      'average_cgd_length' => {},
      'genes_predicted_utr' => {},
      'genome_covered_genes' => {}
    },
    'otherRNA'=> {
      'tRNA_number'=> {},
      'exons_number' => {},
      'total_gene_lenght' => {}
    },
    'function' => {
      'isoform'=>{},
      'PFAM_gene'=>{},
      'INTERPRO_gene'=>{},
      'GO_gene'=>{},
      'KEGG_gene'=>{},
      'MetaCyc_gene'=>{},
      'UniPath_gene'=>{},
      'Reactome_gene'=>{},
      'PFAM_mrna'=>{},
      'INTERPRO_mrna'=>{},
      'GO_mrna'=>{},
      'KEGG_mrna'=>{},
      'MetaCyc_mrna'=>{},
      'UniPath_mrna'=>{},
      'Reactome_mrna'=>{},
      'gene_w_function'=>{},
      'gene_wo_function'=>{},
      'gene_named'=>{},
      'gene_name_duplicated'=>{}
      }
);

#get arg and PARAMETERS
{
  my ($input_file, $output_file);

  # Define script options
  GetOptions(
  'help|h'                => sub { pod2usage( -verbose => 2 ), -exitval => 0 },
  'man'                   => sub { pod2usage( -verbose => 2 )},
  'output-name|o=s'       => \$output_file,
#  'config-name|f=s'       => \$input_file,
  ) or pod2usage(2);

#  pod2usage( "--config-file must be specified" )
#  unless defined $input_file;

  main($input_file, $output_file, %metrics);
}

#MAIN
sub main {
  my ($input_file, $output_file, %metrics) = @_;


  my $ostream     = IO::File->new();
  if(defined($output_file)){
    $ostream->open($output_file, 'w' ) or
    croak(
      sprintf( "Can not open '%s' for reading: %s", $output_file, $! ) );
    }
  else{
    $ostream->fdopen( fileno(STDOUT), 'w' ) or
      croak( sprintf( "Can not open STDOUT for writing: %s", $! ) );
  }


opendir(DIR, $directory) or die $!;

  while (my $file = readdir(DIR)) {

  #  print $file."\n";

        if ($file eq 'config.yaml') {

          $config="yes";
          last;

        }else {
          $config="no";
          #print $config_file;

        }

    }

  if($config eq "yes") {


  #  print $config;
    read_yaml($input_file, $ostream, %metrics);

  } else {

    create_yaml($input_file);

  }

}

sub create_yaml {

  $input_file = "config.yaml";

  my $yaml = IO::File->new();
  if(defined($input_file)){
    $yaml->open($input_file, 'w' ) or
    croak(
      sprintf( "Can not open '%s' for reading: %s", $input_file, $! ) );
  }

  print $yaml "#YAML config files\n";
  print $yaml "#please make sure there is a space between \":\" and your path\n#delete or comment line with # if file not available";
  print $yaml "#write species name eg : Metschnikowia saccharicola\n";
  print $yaml "species: \n";
  print $yaml "#write your specific strain if you have one eg : M_saccharicola\n";
  print $yaml "strain: \n";
  print $yaml "#complete with path and file of your fasta file report eg : \$PATH/genome/M_sacc/fasta_report.txt \n";
  print $yaml "fasta_stats:  \n";
  print $yaml "#complete with path and	file of	your cegma file report eg : \$PATH/Cegma/output.completeness_report\n";
  print $yaml "cegma_stats:  \n";
  print $yaml "#complete with path and	file of	your repeat file report eg : \$PATH/maker/annotations/annotationByType/repeat.txt\n";
  print $yaml "repeat_stats: \n";
  print $yaml "#complete with path and	file of	your gene build file report eg : \$PATH/maker/annotations/annotationByType/MakerManaged/report.txt\n";
  print $yaml "gene_stats: \n";
  print $yaml "#complete with path and	file of	your functional annotation file report eg : \$PATH/functional_annotation/codingGeneFeatures/codingGeneFeatures_withAnnotation.gff/report.txt\n";
  print $yaml "function_stats: \n";
  print $yaml "#complete with path and	file of	your RNA file report eg : \$PATH/maker/annotations/annotationByType/MakerManaged/report.txt\n";
  print $yaml "otherRNA_stats: \n";

  print
  "  -------------------------------------------------------------------\n
  \t\tconfig.yaml file has been created\n
  \t\tplease fill it with correct paths\n
  \t\tthen rerun the script to generate the report. \n
  -------------------------------------------------------------------\n";

}



sub read_yaml {
  #get files from YAML config file
  my ($input_file, $ostream, %metrics) = @_;

  $input_file = 'config.yaml';

  $config_file = LoadFile($input_file);

  #check if files exist and run corresponding functions
  if(exists $config_file->{'species'}){
    $species = $config_file->{'species'};
  }else{
    print "There is no species provided in the YamL config file\n";
  }
  if(exists $config_file->{'strain'}){
  $strain = $config_file->{'strain'};
  }else{
    print "There is no strain provided in the YamL config file\n";
  }
  if(exists $config_file->{'fasta_stats'}){
    $fasta_file = $config_file->{'fasta_stats'};
    read_fasta($fasta_file, %metrics);
  }else{
    print "There is no fasta file report provided in the YamL config file\n";
  }
  if(exists $config_file->{'repeat_stats'}){
    $repeat_file = $config_file->{'repeat_stats'};
    read_repeats($repeat_file, %metrics);
  }else{
    print "There is no repeat file report provided in the YamL config file\n";
  }
  if(exists $config_file->{'gene_stats'}){
    $gene_file = $config_file->{'gene_stats'};
    read_gene($gene_file, %metrics);
  }else{
    print "There is no gene file report provided in the YamL config file\n";
  }
  if(exists $config_file->{'otherRNA_stats'}){
    $otherRNA_file = $config_file->{'otherRNA_stats'};
    read_otherRNA($otherRNA_file, %metrics);
  }else{
    print "There is no otherRNA file report provided in the YamL config file\n";
  }
  if(exists $config_file->{'cegma_stats'}){
    $cegma_file = $config_file->{'cegma_stats'};
    read_cegma($cegma_file, %metrics);
  }else{
    print "There is no cegma file report provided in the YamL config file\n";
  }
  if(exists $config_file->{'function_stats'}){
    $function_file = $config_file->{'function_stats'};
    read_function($function_file, %metrics);
  }else{
    print "There is no functional file report provided in the YamL config file\n";
  }
  print_hash($ostream, $config_file, %metrics);
}

#FUNCTIONS
#get info from fasta file
sub read_fasta {
  my ($fasta_file, %metrics) = @_;

  my $fasta = IO::File->new();
  if ( defined $fasta_file ) {
    $fasta->open($fasta_file, 'r') or croak (sprintf("Can not open '%s' for reading: %s", $fasta_file, $!));
  }

  #get information
  while(my $line =<$fasta>) {
    chomp $line;
   if ($line=~/(\d+)\ssequences/){
     $metrics{'fasta'}->{'seq_number'}=$1;
    }
    if ($line=~/(\d+)\snucleotides\D+(\d+)/){
      $metrics{'fasta'}->{'nuc_number'}=$1;
      $metrics{'fasta'}->{'N_number'}=$2;
    }
    if ($line=~/(\d+)\sN-regions/){
       $metrics{'fasta'}->{'N_region_number'}=$1;
     }
    if ($line=~/(\d+.\d+%)/){
      $metrics{'fasta'}->{'GC_content'}=$1;
    }
    if ($line=~/N50\D+(\d+)$/){
      $metrics{'fasta'}->{'N50'}=$1;
    }
  }
}
#get info from cegma file
sub read_cegma {
  my ($cegma_file, $ostream) = @_;

  my $cegma = IO::File->new();
  if ( defined $cegma_file ) {
    $cegma->open($cegma_file, 'r') or croak (sprintf("Can not open '%s' for reading: %s", $cegma_file, $!));
  }
  #get information
  while(my $line =<$cegma>) {
    chomp $line;
    if ($line=~/Complete[\s\t]+(\d+)[\s\t]+(\d*\.*\d*)[\s\t]+\-[\s\t]+(\d+)[\s\t]+(\d*\.*\d*)[\s\t]+(\d*\.*\d*)/){
      $metrics{'cegma'}->{'prot_complete_number'}=$1;
      $metrics{'cegma'}->{'percent_completeness'}=$2;
      $metrics{'cegma'}->{'complete_total'}=$3;
      $metrics{'cegma'}->{'complete_average'}=$4;
      $metrics{'cegma'}->{'percent_complete_orthologs'}=$5;
    }
    if ($line=~/Partial[\s\t]+(\d+)[\s\t]+(\d*\.*\d*)[\s\t]+\-[\s\t]+(\d+)[\s\t]+(\d*\.*\d*)[\s\t]+(\d*\.*\d*)/){
      $metrics{'cegma'}->{'prot_partial_number'}=$1;
      $metrics{'cegma'}->{'percent_partial'}=$2;
      $metrics{'cegma'}->{'partial_total'}=$3;
      $metrics{'cegma'}->{'partial_average'}=$4;
      $metrics{'cegma'}->{'percent_partial_orthologs'}=$5;
    }
  }
}
#get info from repeat file
sub read_repeats {
  my ($repeat_file, $ostream) = @_;

  my $repeat = IO::File->new();
  if ( defined $repeat_file ) {
    $repeat->open($repeat_file, 'r') or croak (sprintf("Can not open '%s' for reading: %s", $repeat_file, $!));
  }
  #get information
  while(my $line =<$repeat>) {
    chomp $line;
    if ($line=~/Total\t(\d+)\t(\d*\.*\d*)\t(\d*\.*\d*)\t(\d*\.*\d*)/){
      $metrics{'repeats'}->{'number'}=$1;
      $metrics{'repeats'}->{'total_size'}=$2;
      $metrics{'repeats'}->{'mean_size'}=$3;
      $metrics{'repeats'}->{'percent_genome'}=$4;
    }
  }
}
#get info from gene file and info from other RNA (for now only tRNA)
sub read_gene {

  my $gene_boolean;
  my $UTR_nb;
  my $mRNA_nb;
  my $fraction_utr_mRNA;

  my ($gene_file, $ostream) = @_;

  my $gene = IO::File->new();
  if ( defined $gene_file ) {
    $gene->open($gene_file, 'r') or croak (sprintf("Can not open '%s' for reading: %s", $gene_file, $!));
  }
  #get information
  while(my $line =<$gene>) {
    chomp $line;

    #OTHER RNA
    if ($line=~/Information about Non_Coding_Gene/){
      $gene_boolean = 'false';
    }
    #gene_build
    if ($line=~/Information about Coding_Gene/){
      $gene_boolean = 'true';
    }

    if ($line=~/Number of genes[\s\t]+(\d+)/) {

      if ($gene_boolean eq 'true') {
        $metrics{'gene_build'}->{'prot_coding_gene_nb'}=$1;
      }
    }

    if ($line=~/mean cds length[\s\t]+(\d+)/) {
      $metrics{'gene_build'}->{'average_cgd_length'}=$1;
    }


     if ($line=~/Number of mrnas with at least one utr[\s\t]+(\d+)/) {
       $UTR_nb = $1;
     }

     if ($line=~/Number of mrnas[\s\t]+(\d+)/) {
       $mRNA_nb = $1;
     }

    if ($line=~/% of genome covered by gene[\s\t]+(\d*\.*\d*)/) {
      if ($gene_boolean eq 'true') {
        $metrics{'gene_build'}->{'genome_covered_genes'}=$1;
      }
    }
  }
  $fraction_utr_mRNA = sprintf("%.2f", ($UTR_nb*100)/$mRNA_nb);
 #$fraction_utr_mRNA = ($UTR_nb*100)/$mRNA_nb;
 $metrics{'gene_build'}->{'genes_predicted_utr'}=$fraction_utr_mRNA;

}

sub read_otherRNA {

  my $gene_boolean;

  my ($otherRNA_file, $ostream) = @_;

  my $otherRNA = IO::File->new();
  if ( defined $otherRNA_file ) {
    $otherRNA->open($otherRNA_file, 'r') or croak (sprintf("Can not open '%s' for reading: %s", $otherRNA_file, $!));
  }
  #get information
  while(my $line =<$otherRNA>) {
    chomp $line;

    #OTHER RNA
    if ($line=~/Information about Non_Coding_Gene/){
      $gene_boolean = 'false';
    }
    #gene_build
    if ($line=~/Information about Coding_Gene/){
      $gene_boolean = 'true';
    }

    if ($line=~/Number of genes[\s\t]+(\d+)/) {

      if ($gene_boolean eq 'false') {
        $metrics{'otherRNA'}->{'tRNA_number'}=$1;
      }
    }

    if ($line=~/Number of exons[\s\t]+(\d+)/) {
      if($gene_boolean eq 'false'){
      $metrics{'otherRNA'}->{'exons_number'}=$1;
      }
    }
    if ($line=~/Total gene length[\s\t]+(\d+)/) {
      if($gene_boolean eq 'false'){
      $metrics{'otherRNA'}->{'total_gene_lenght'}=$1;
      }
    }
  }
}
#get info from functional annotation file
sub read_function {
  my ($function_file, $ostream) = @_;

  my $function = IO::File->new();
  if ( defined $function_file ) {
    $function->open($function_file, 'r') or croak (sprintf("Can not open '%s' for reading: %s", $function_file, $!));
  }
  #get information
  while(my $line =<$function>) {
    chomp $line;

    if ($line=~/PFAM[\s\t]*\|/){
      my @pfam = split(/[\s\t]*\|/, $line);
      $metrics{'function'}->{'PFAM_gene'}=$pfam[5];
      $metrics{'function'}->{'PFAM_mrna'}=$pfam[4];
    }
    if ($line=~/InterPro[\s\t]*\|/){
      my @InterPro = split(/[\s\t]*\|/, $line);
      $metrics{'function'}->{'INTERPRO_gene'}=$InterPro[5];
      $metrics{'function'}->{'INTERPRO_mrna'}=$InterPro[4];
    }

    if ($line=~/GO[\s\t]*\|/){
      my @GO = split(/[\s\t]*\|/, $line);
      $metrics{'function'}->{'GO_gene'}=$GO[5];
      $metrics{'function'}->{'GO_mrna'}=$GO[4];
    }

    if ($line=~/KEGG[\s\t]*\|/){
      my @KEGG = split(/[\s\t]*\|/, $line);
      $metrics{'function'}->{'KEGG_gene'}=$KEGG[5];
      $metrics{'function'}->{'KEGG_mrna'}=$KEGG[4];
    }

    if ($line=~/MetaCyc[\s\t]*\|/){
      my @MetaCyc = split(/[\s\t]*\|/, $line);
      $metrics{'function'}->{'MetaCyc_gene'}=$MetaCyc[5];
      $metrics{'function'}->{'MetaCyc_mrna'}=$MetaCyc[4];
    }

    if ($line=~/UniPathwa\|/){
      #print $line."\n";
      my @UniPath = split(/[\s\t]*\|/, $line);
    #  print $UniPath[1]."\n";
      $metrics{'function'}->{'UniPath_gene'}=$UniPath[5];
      $metrics{'function'}->{'UniPath_mrna'}=$UniPath[4];
    #  print Dumper \@UniPath;
    }

    if ($line=~/Reactome[\s\t]*\|/){
      my @Reactome = split(/[\s\t]*\|/, $line);
      $metrics{'function'}->{'Reactome_gene'}=$Reactome[5];
      $metrics{'function'}->{'Reactome_mrna'}=$Reactome[4];
    }


    if($line=~/(\d+) mRNA isoforms/){
      $metrics{'function'}->{'isoforms'}=$1;
    }
    if ($line=~/nb gene with\s+\D+(\d+)/) {
        $metrics{'function'}->{'gene_w_function'}=$1;
    }

    if ($line=~/nb gene without\s+\D+(\d+)/) {
        $metrics{'function'}->{'gene_wo_function'}=$1;
    }
    if ($line=~/(\d+)\D+names\D+(\d+)\D+duplicated/){
      $metrics{'function'}->{'gene_named'}=$1;
      $metrics{'function'}->{'gene_name_duplicated'}=$2;
    }
  }
}
#function to print everything
sub print_hash {

  my ($ostream, $config_file, %metrics) = @_;

  if(exists $config_file->{'species'}){
  print $ostream "\nSPECIES\t".$config_file->{'species'}."\n\n";
  } else{
  print $ostream "\nThere is no species information available\n";
  }
  # print fasta info
  if(exists $config_file->{'fasta_stats'}){
    foreach my $fasta_info ('fasta') {
      print $ostream "-------------------------------\nGENOME ASSEMBLY\n-------------------------------\n";
      print $ostream "Strains\tSequence number(over 1kb)\tNucleotide number\tN number\tN-regions number\tGC-content\tN50\n";
      print $ostream $config_file->{'strain'}."\t";
      print $ostream $metrics{$fasta_info}->{'seq_number'}."\t";
      print $ostream $metrics{$fasta_info}->{'nuc_number'}."\t";
      print $ostream $metrics{$fasta_info}->{'N_number'}."\t";
      print $ostream $metrics{$fasta_info}->{'N_region_number'}."\t";
      print $ostream $metrics{$fasta_info}->{'GC_content'}."\t";
      print $ostream $metrics{$fasta_info}->{'N50'}."\n";
      print $ostream "\n\n";
    }
  }else{
    print $ostream "There is no Genome assembly information available\n";
  }

  #print cegma info
  if(exists $config_file->{'cegma_stats'}){
    foreach my $cegma_info ('cegma') {
      print $ostream "-------------------------------\nCEGMA\n-------------------------------\n";
      print $ostream "Strain\tProtein Number\tPercentage Completeness\tTotal\tAverage\tPercentage Orthologs\n";
      print $ostream $config_file->{'strain'}."\t";
      print $ostream $metrics{$cegma_info}->{'prot_complete_number'}."\t";
      print $ostream $metrics{$cegma_info}->{'percent_completeness'}."%\t";
      print $ostream $metrics{$cegma_info}->{'complete_total'}."\t";
      print $ostream $metrics{$cegma_info}->{'complete_average'}."\t";
      print $ostream $metrics{$cegma_info}->{'percent_complete_orthologs'}."%\n";
      print $ostream "\t";
      print $ostream $metrics{$cegma_info}->{'prot_partial_number'}."\t";
      print $ostream $metrics{$cegma_info}->{'percent_partial'}."%\t";
      print $ostream $metrics{$cegma_info}->{'partial_total'}."\t";
      print $ostream $metrics{$cegma_info}->{'partial_average'}."\t";
      print $ostream $metrics{$cegma_info}->{'percent_partial_orthologs'}."%\n";
      print $ostream "\n\n";
    }
  }else{
    print $ostream "There is no cegma information available\n";
  }
  #print repeat info
  if(exists $config_file->{'repeat_stats'}){
    foreach my $repeat_info ('repeats') {
      print $ostream "-------------------------------\nREPEAT MASKING\n-------------------------------\n";
      print $ostream "Strain\tNumber\tTotal size(kb)\tMean size(bp)\t% genome\n";
      print $ostream $config_file->{'strain'}."\t";
      print $ostream $metrics{$repeat_info}->{'number'}."\t";
      print $ostream $metrics{$repeat_info}->{'total_size'}."\t";
      print $ostream $metrics{$repeat_info}->{'mean_size'}."\t";
      print $ostream $metrics{$repeat_info}->{'percent_genome'}."%\n";
      print $ostream "\n\n";
    }
  }else{
    print $ostream "There is no repeat information available\n";
  }
  #print gene info
  if(exists $config_file->{'gene_stats'}){
    foreach my $gene_info('gene_build') {
      print $ostream "-------------------------------\nGENE BUILD\n-------------------------------\n";
      print $ostream "Strain\tNumber of protein-coding genes\tAverage CDS lenght\tFraction of genes with predicted UTR\tFraction of the genome covered by genes\n";
      print $ostream $config_file->{'strain'}."\t";
      print $ostream $metrics{$gene_info}->{'prot_coding_gene_nb'}."\t";
      print $ostream $metrics{$gene_info}->{'average_cgd_length'}."\t";
      print $ostream $metrics{$gene_info}->{'genes_predicted_utr'}."%\t";
      print $ostream $metrics{$gene_info}->{'genome_covered_genes'}."%\n";
      print $ostream "\n\n";
    }
  }else{
    print $ostream "There is no gene information available\n";
  }
  #print function info
  if(exists $config_file->{'function_stats'}){
    foreach my $function_info('function') {
      print $ostream "-------------------------------\nFUNCTIONAL ANNOTATION\n-------------------------------\n";
      print $ostream "Strain\tPFAM\tInterpro\tGO\tKEGG\tMetaCyc\tUnipathway\tReactome\n";
      print $ostream $config_file->{'strain'}."\t";
      print $ostream $metrics{$function_info}->{'PFAM_gene'}."\t";
      print $ostream $metrics{$function_info}->{'INTERPRO_gene'}."\t";
      print $ostream $metrics{$function_info}->{'GO_gene'}."\t";
      print $ostream $metrics{$function_info}->{'KEGG_gene'}."\t";
      print $ostream $metrics{$function_info}->{'MetaCyc_gene'}."\t";
      print $ostream $metrics{$function_info}->{'UniPath_gene'}."\t";
      print $ostream $metrics{$function_info}->{'Reactome_gene'}."\n";

      if ($metrics{$function_info}->{'isoforms'} != '0') {
        print $ostream "\t";
        print $ostream $metrics{$function_info}->{'PFAM_mrna'}."\t";
        print $ostream $metrics{$function_info}->{'INTERPRO_mrna'}."\t";
        print $ostream $metrics{$function_info}->{'GO_mrna'}."\t";
        print $ostream $metrics{$function_info}->{'KEGG_mrna'}."\t";
        print $ostream $metrics{$function_info}->{'MetaCyc_mrna'}."\t";
        print $ostream $metrics{$function_info}->{'UniPath_mrna'}."\t";
        print $ostream $metrics{$function_info}->{'Reactome_mrna'}."\n";
      }
      print $ostream "\n\n";
      print $ostream "Strain\tGene with functional annotation\tGene without functional annotation\n";
      print $ostream $config_file->{'strain'}."\t";
      print $ostream $metrics{$function_info}->{'gene_w_function'}."\t";
      print $ostream $metrics{$function_info}->{'gene_wo_function'}."\n";

      print $ostream "\n\n";
      print $ostream "-------------------------------\nGENE NAME INFERENCE\n-------------------------------\n";
      print $ostream "Strain\tGene named\tNumber of duplicated name among gene named\n";
      print $ostream $config_file->{'strain'}."\t";
      print $ostream $metrics{$function_info}->{'gene_named'}."\t";
      print $ostream $metrics{$function_info}->{'gene_name_duplicated'}."\n";
      print $ostream "\n\n";
    }
  }else{
    print $ostream "There is no functional information available\n";
  }
  #print other RNA
  if(exists $config_file->{'otherRNA_stats'}){
    foreach my $otherRNA_info('otherRNA') {
      print $ostream "-------------------------------\ntRNA\n-------------------------------\n";
      print $ostream "Strain\tNumber of tRNA\tNumber of exons\tTotal gene lenght\n";
      print $ostream $config_file->{'strain'}."\t";
      print $ostream $metrics{$otherRNA_info}->{'tRNA_number'}."\t";
      print $ostream $metrics{$otherRNA_info}->{'exons_number'}."\t";
      print $ostream $metrics{$otherRNA_info}->{'total_gene_lenght'}."\n";
      print $ostream "\n\n";
    }
  }else{
    print $ostream "There is no other RNA information available\n";
  }
}

#print Dumper(%metrics);

1;

__END__;

=head1 NAME

create a pre tabulated report for the annotation final report

=head1 AUTHOR

Lucile SOLER NBIS 29/02/2016

=head1 SYNOPSIS

create pretabulated report

report.pl --help

report.pl --input-name

report.pl --output-name

=head1 OPTIONS

=over

=item B<--help>

Display a brief usage message.

=item B<--man>

Display the manual page.

=item B<--input-file|-f>

Give the config file where all paths of all different outputs are required for report

=item B<--output-file|-o>

Give the output file

=back

=head1 DESCRIPTION

create a pre tabulated report for the annotation final report

=cut
