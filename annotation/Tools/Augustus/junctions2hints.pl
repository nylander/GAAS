#!/usr/bin/perl
# A script to convert a junction.bed output from Tophat to a gff format to be used as an intronhints file in Augustus
# If you have several junctions.bed files, concatenate these first
# usage: perl junctions2hints.pl --infile junctions.bed


use warnings;
use strict;
use Getopt::Long;

my $usage = qq{
Usage: perl junctions2hints.pl --infile junctions.bed
If you have several junctions.bed files, concatenate these first.
A file called intronhints.gff will be created with all your intron-hints, including multiplicities.

 Getting help:
    [--help]

  Input file
    [--infile]
    Name of junctions.bed file
  Output file
    [--outfile]
    Name of the hint file to write
};

my $help;
my $infile;
my $outfile;
my %junctions=();

GetOptions(
    "help" => \$help,
    "outfile=s" => \$outfile,
    "infile=s" => \$infile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}
open INFILE, $infile or die "$usage";
open HINTS,">$outfile";

while (<INFILE>) {
  chomp;
  unless ($_ =~ /^track/){
    my @bed_line=split(/\t/, $_);
    my ($startblock,$endblock)=split(/\,/, $bed_line[10]);
    $bed_line[1]=$bed_line[1]+$startblock+1;
    $bed_line[2]=$bed_line[2]-$endblock;
    my $key = join (':',$bed_line[0],$bed_line[1],$bed_line[2]);
    $junctions{$key} +=$bed_line[4];
  }
}
close INFILE;
 
foreach my $key (sort keys %junctions){
  my ($scaffold, $start, $stop) = split (/\:/, $key);
  print HINTS "$scaffold\ttophat\tintron\t$start\t$stop\t0\t\.\t\.\tmult=$junctions{$key}\;src=E\n";
}

close HINTS;
