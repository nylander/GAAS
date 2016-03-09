#!/usr/bin/perl

## BILS 2015 (www.bils.se)
## jacques.dainat@bils.se

use strict;
use Getopt::Long;
use Bio::Tools::GFF;
use Pod::Usage;

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--gff filename]
		The name of the cegma gff file to convert. 
		
  Ouput:    
    [--outfile filename]
        The name of the output file (A bed file).

  At least the input cegma gff file is mandatory:
  Usage: script.pl --gff infile.gff [--outfile outfile.bed]


  This script allows to convert a gff file from cegma output to a bed file.
};

my $outfile = undef;
my $gff = undef;
my $attributes = undef ;
my $help;

GetOptions(
    "help" => \$help,
    "gff=s" => \$gff,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if ( ! (defined($gff)) ){
    print $usage;
    exit(0);
}

## Manage output file
my $fh;
if ($outfile) {
open($fh, '>', $outfile) or die "Could not open file '$outfile' $!";
}
else{ 
  $fh=\*STDOUT;
}

### Read gff input file 
open(my $fhIn, '<:encoding(UTF-8)', $gff)
  or die "Could not open file '$gff' $!";
 
## Read line by line and save in hash with geneName as key. Allows to merge feature belonging to the same gene.
my %hash;
while (my $row = <$fhIn>) {
  chomp $row;
  my @tab = split( /\s/,$row);
  if (lc($tab[2]) eq "exon"){
    if(exists($hash{$tab[8]})){
        push (@{$hash{$tab[8]}}, $row)
    }
    else{$hash{$tab[8]}=[$row]}  
  }
}
$fhIn->close();


my $nbGene=keys %hash;
print "$nbGene genes read.\n";

# foreach gene we have a list of feature
foreach my $key (keys %hash){
  my @tabValues=@{$hash{$key}};
  my $nbValue=0;
  my @listSizeStart;
  #foreach fetaure of a gene
  foreach my $value (@tabValues){
    #cut the gff feature
    my ($chr,$tool,$type,$start,$stop,$score,$dir,$col8,$name) = split(/\s/,$value);
    my $size=$stop-$start;
    #In the case where we have only one feature
    if ($#tabValues == 0){ 
      print $fh "$chr\t$start\t$stop\t$name\t$score\t$dir\t$start\t$stop\t0\t1\t$size,\t0\n";
    }
    #In the case where we have several features we save information
    elsif ($nbValue!=$#tabValues){
      push(@listSizeStart,[$size,$start]);
      $nbValue++;
    }
    #In the case where we have several features we manage information kept and manage them to print the good final bed feature.
    else{
      push(@listSizeStart,[$size,$start]); #save information needed of the last feature of the gene
      my (@listSorted)=sortByPos(@listSizeStart); #sort the infromation
      my $final_sizeList; my $final_startList;
      my $cpt=0;my $originStart;
      #foreach information save we will create correct strings for the output
      foreach my $tabDuoRef (@listSorted){
        my ($res_size,$res_start)=@{$tabDuoRef};
        if($cpt==0){ $originStart=$res_start; $cpt++;}
          my $start_corrected=$res_start-$originStart;
          $final_sizeList.="$res_size,";
          $final_startList.="$start_corrected,";
      }
      #print result in case where gene has several features
      print $fh "$chr\t$start\t$stop\t$name\t$score\t$dir\t$start\t$stop\t0\t1\t$final_sizeList\t$final_startList\n";
    }

  }
}
$fh->close();




#This function sort a tab of tab according to the second value of the sub-tab ...
sub sortByPos {
  my (@featureList) = @_;

  my @featureListSorted =  ( sort ({ $a->[1] <=> $b->[1] } @featureList));
  return @featureListSorted;
} 
