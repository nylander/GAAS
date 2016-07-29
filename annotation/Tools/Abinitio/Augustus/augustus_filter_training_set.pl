#!/usr/bin/perl

use strict;
use Getopt::Long;
use Scalar::Util qw(openhandle);
use Bio::
use Time::Piece;
use Time::Seconds;

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--gff filename]
		The name of the file to read.
		
	[--genome filename]
		The genome sequence to read.
		
  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile = undef;
my $gff = undef;
my $genome = undef;
my $help;

GetOptions(
    "help" => \$help,
    "gff=s" => \$gff,
	"genome=s" => \$genome,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if ($outfile) {
    open(STDOUT, ">$outfile") or die("Cannot open $outfile");
}

# Create protein dataset

run_blast(gff,genome);


#my $gffio = Bio::Tools::GFF->new(-file => $gff, -gff_version => 3);







sub run_blast {
	
	my $gff = shift;
	my $genome = shift;
	
	print $gff ;
	
}


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


