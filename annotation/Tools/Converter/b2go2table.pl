#!/usr/bin/perl

use strict;
use Getopt::Long;
use Bio::Tools::GFF;


my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:

        [--b2go filename]
                The name of the Blast2Go annotation file to read

  Ouput:
    [--outfile filename]
        The name of the output file.
};

my $outfile = undef;
my $b2go = undef;
my $help;

GetOptions(
    "help" => \$help,
        "b2go=s" => \$b2go,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if ($outfile) {
    open(STDOUT, ">$outfile") or die("Cannot open $outfile");
}


my %lookup = {};

open (my $IN, '<', $b2go) or die "FATAL: Can't open file: $b2go for reading.\n$!\n";

while (<$IN>) {
        chomp;
        my $line = $_;
        my ($id,$go,$name) = split("\t", $line);
        if (exists $lookup{$id}) {
                $lookup{$id} .= "," . $go ;
        } else {
                $lookup{$id} = $go ;
        }

}

foreach my $transcript (keys %lookup) {

        my $go_terms = $lookup{$transcript};

        print $transcript . "\t" . $go_terms . "\n";
}

close ($IN);
	



