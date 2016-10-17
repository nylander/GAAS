#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Scalar::Util qw(openhandle);
use Time::Piece;
use Time::Seconds;
use Cwd;

my $dir = getcwd;

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]
  	
  Ouput:    
    [--outfile filename]
        The name of the masked genome file. By default, the name will genome.rm.fa

   Usage:
	Must be executed in the folder from which Maker was run and will find the maker output
	on its own and create a concatenated masked assembly. 

};

my $outfile = "genome.rm.fa";
my $maker_dir = undef;
my $datastore = undef;
my $quiet;
my $help;

GetOptions(
    "help" => \$help,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if (-f $outfile) {
	die "The outfile $outfile already exists, exiting\n";
} 

# Find the datastore index

my $maker_dir = undef;

opendir(DIR, $dir) or die "couldn't open $dir: $!\n";
my @files = readdir DIR;
closedir DIR;

if (my ($matched) = grep $_ =~ /^.*\.maker\.output$/ , @files) {
    $maker_dir = $dir . "/" . $matched ;
    my $base = $matched;
    $base =~ s/\.maker\.output//g ;
    $datastore = $matched . "/" . $base . "_master_datastore_index.log" ;
}

die "There seems to be no maker output directory here, exiting...\n" unless defined $maker_dir;

if (-f $dir . "/" . $datastore ) {
        print "Found datastore, merging masked contigs...\n";
} else {
        die "Could not find datastore index ($datastore), exiting...\n";
}



# This is one way to open a file...
open (my $IN, '<', $datastore) or die "FATAL: Can't open file: $datastore for reading.\n$!\n";

# Streaming the file, line by line
while (<$IN>) {
	chomp; # Trims the line, removes the line breaks
	my $line = $_; # store the line in a variable
	next unless ($line =~ /^.*FINISHED$/) ; # We only want finished contig annotations...
	
	my ($contig,$location,$status) = split("\t",$line);
	
	# If the contig includes a dot character, the output file will include a percent character...
	if ($contig =~ /^.*\..*$/) {
		$contig =~ s/\./\%2E/g ;
	}

	my $void = $location . "theVoid." . $contig ;
	
	my $repeat_contig = $maker_dir . "/" . $void . "/" . "query.masked.fasta" ;
	
	if (-f $repeat_contig ) {
		system("cat $repeat_contig >> $outfile");
	}
	
}
# We should close the file to make sure that the transaction finishes cleanly.
close ($IN);


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


