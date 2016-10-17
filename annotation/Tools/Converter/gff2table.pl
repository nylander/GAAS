#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Scalar::Util qw(openhandle);
use Bio::Tools::GFF;
use Time::Piece;
use Time::Seconds;

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--gff filename]
		The name of the file to read. 
  Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile = undef;
my $gff = undef;
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

if ($outfile) {
    open(STDOUT, ">$outfile") or die("Cannot open $outfile");
}

open(GFF, "<$gff") || die("Can't open $gff.");

my $gffio = Bio::Tools::GFF->new(-file => $gff, -gff_version => 3);
my $gffout = Bio::Tools::GFF->new(-gff_version => 3);


while( my $feature = $gffio->next_feature()) {
	
	if ( $feature->primary_tag =~ /gene/ or $feature->primary_tag =~ /mRNA/ ) {
		
		my $description = "" ;
		my $dbxref = "" ;
		
		if ($feature->primary_tag eq "mRNA") {
			
			if ($feature->has_tag('description') ) {
				my @values = $feature->get_tag_values('description');
				$description = join(",",@values) ;
			}
			
			if ($feature->has_tag('Dbxref') ) {
				my @values = $feature->get_tag_values('Dbxref') ;
				$dbxref = join(",",@values) ;
			}
			
		}	
			
		my @id_values = $feature->get_tag_values('ID');
		my $id = shift @id_values;
		
		
		print $feature->primary_tag . "\t" . $id . "\t" . $feature->seq_id . "\t" . $feature->start . "\t" . $feature->end . "\t" . $feature->strand . "\t" . $description . "\t" . $dbxref . "\n" ;
	
	}
	
}

$gffio->close();


# --------------

sub msg {
  my $t = localtime;
  my $line = "[".$t->hms."] @_\n";
}

sub runcmd {
  msg("Running:", @_);
  system(@_)==0 or err("Could not run command:", @_);
}

sub err {
  msg(@_);
  exit(2);
}


