#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Scalar::Util qw(openhandle);
use Time::Piece;
use Time::Seconds;

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--infile filename]
		The name of the file to read. 
	[--type biotype]
		Indicate whether input sequences are protein (p) or nucleotide (n). Default is n, nucleotide.
  Ouput:    
    [--outdir foldername]
        The name of the output folder.
};

my $outdir = undef;
my $infile = undef;
my $type = "n";
my $quiet = undef;
my $help;

GetOptions(
    "help" => \$help,
	"type=s" => \$type,
    "infile=s" => \$infile,
    "outdir=s" => \$outdir);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

defined $outdir or die "Must specify a directory (--outdir)\n";
defined $infile or die "Must specify an input fasta file (--infile)\n";

# Create outdir if it does not consist

if (-d $outdir) {
     msg("Outdir found!");
} elsif (-f $outdir . "/" . $infile . ".gff3") {
	msg("This folder already contains an InterProScan analysis. Aborting!") and die;
} else {
	runcmd("mkdir $outdir");
	msg("Outdir created!");
}

msg("Starting InterProScan");

# Check whether sequence is nucleotide or protein

# Run InterProScan

runcmd("/sw/bioinfo/interproscan-5.2-45.0/interproscan.sh -appl PfamA-27.0,TIGRFAM-13.0,ProDom-2006.1 -d $outdir -i $infile -f GFF3,XML -iprlookup -ms 40 -T /tmp -t $type");

# --------------

sub msg {
  my $t = localtime;
  my $line = "[".$t->hms."] @_\n";
  print LOG $line if openhandle(\*LOG);
  print STDERR $line unless $quiet;
}

sub runcmd {
  msg("Running:", @_);
  system(@_)==0 or err("Could not run command:", @_);
}

sub err {
  $quiet=0;
  msg(@_);
  exit(2);
}


