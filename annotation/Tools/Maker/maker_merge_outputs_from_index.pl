#!/usr/bin/perl

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
    [--outdir filename]
        The name of the output directory. By default the file is named annotations.gff

   Usage:
	Must be executed in the folder from which Maker was run and will find the maker output
	on its own and create a concatenated annotation file. 

};
my $outdir = "annotations";
my $protein_file = "annotations.proteins.fa";
my $annotations_file = "annotations.gff";
my $maker_dir = undef;
my $datastore = undef;
my $quiet;
my $help;

GetOptions(
    "help" => \$help,
    "outdir=s" => \$outdir);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if (-d "$outdir") {
	die "The outdirectory $outdir already exists, exiting\n";
} 
else{
	mkdir $outdir;
}

# Find the datastore index

my $maker_dir = undef;

opendir(DIR, $dir) or die "couldn't open $dir: $!\n";
my @files = readdir DIR;
closedir DIR;

my (@matchedDir) = grep $_ =~ /^.*\.maker\.output$/ , @files ;
my $nbDir=$#matchedDir+1;
print "We found $nbDir maker output directorie(s):\n";
foreach my $j (@matchedDir){
print "\t+$j\n";}
if ($nbDir == 0){die "There seems to be no maker output directory here, exiting...\n";}
elsif ($nbDir > 1 ){print "Results will be merged together !\n";}


foreach my $matched (@matchedDir){

	$maker_dir = $dir . "/" . $matched ;
	my $base = $matched;
	$base =~ s/\.maker\.output//g ;
	$datastore = $matched . "/" . $base . "_master_datastore_index.log" ;

	if (-f $dir . "/" . $datastore ) {
        	print "Found datastore in $matched, merging annotations now...\n";
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

		####
		#NOTE: Not all special character will be in URI format. As example the underscore.

		# If the contig includes a dot character, the output file will include a percent character...
		if ($contig =~ /^.*\..*$/) {
			$contig =~ s/\./\%2E/g ;
		}
	
		# If the contig includes a pipe character, the output file will include a percent character...
   		if ($contig =~ m/\|/) {
    		$contig =~ s/\|/\%7C/g ;
   		}

		my $gff_file = $maker_dir . "/" . $location . $contig . ".gff" ;
		my $aa_file = $maker_dir . "/" . $location . $contig . ".maker.proteins.fasta" ;

		if (-f $gff_file ) {
			system("cat $gff_file >> $outdir/$annotations_file");
		}
	
		if (-f $aa_file ) {
			system("cat $aa_file >> $outdir/$protein_file");
		}
	}
	# We should close the file to make sure that the transaction finishes cleanly.
	close ($IN);
}

# Now manage to split file by kind of data
print "Now split file by data type...\n";
my $splitedData_dir= "$outdir/annotationByType";
mkdir $splitedData_dir;

#call split script
my $SplitScript="/projects/scripts/gmod/split_gff_by_source.pl";
system("$SplitScript","--input","$outdir/$annotations_file","-d","$splitedData_dir");


#convert the gff in gtf
if (-f "${splitedData_dir}/maker.gff"){
	print "Converting Maker file to GTF...\n";
	my $gffreadPath="/sw/bioinfo/cufflinks/cufflinks-2.1.1/gffread";
	system("$gffreadPath","-o","$splitedData_dir/maker.gtf","-T","-F","$splitedData_dir/maker.gff");
}
else{print "No gff file to convert\n";}
print "All done!\n";

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


