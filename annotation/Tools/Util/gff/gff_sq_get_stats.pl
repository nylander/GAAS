#!/usr/bin/perl

use strict;
use Getopt::Long;
use Scalar::Util qw(openhandle);
use Time::Piece;
use Time::Seconds;
use List::MoreUtils 'any';

use Bio::Tools::GFF;


my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--infile filename]
		The name of the GFF/GTF file to read. 
	
};

my $outfile = undef;
my $infile = undef;
my $help;

GetOptions(
    "help" => \$help,
    "infile=s" => \$infile,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

if ($outfile) {
    open(STDOUT, ">$outfile") or die("Cannot open $outfile");
}

# Here we store all the numbers
my %data_store = ( 'exons' => 0 , 'mRNAs' => 0 , 'genes' => 0 , 'coding_nt' => 0  );

# We want to find out what sort of file we are presented with:
my $gff_file;
my $is_gtf = 0;

if ($infile =~ /\.gff/) {
	$gff_file = Bio::Tools::GFF->new(-file =>$infile , -gff_version => 3);
} elsif ($infile =~ /\.gtf/) {
	$is_gtf = 1;
	$gff_file = Bio::Tools::GFF->new(-file => $infile , -gff_version => 2.5 );
} else {
		die "Sorry, no clue what file format this is (need extension .gtf or .gff)\n";
}

my @transcript_ids = () ;
my @gene_ids = ();

while(my $feature = $gff_file->next_feature()) {
    
	if ($feature->primary_tag eq "exon" ) {
		$data_store{'exons'} += 1 ;
	} elsif ($feature->primary_tag =~ /[mstsno]RNA/ or $feature->primary_tag eq "transcript") {
		$data_store{'mRNAs'} += 1 ;
	} elsif($feature->primary_tag eq "gene") {
		$data_store{'genes'} += 1;
	} elsif ($is_gtf == 1) {
		
		my @tvalues = $feature->get_tag_values('transcript_id');
		my $transcript_id = shift @tvalues;
		
		unless (my ($matched) = grep $_ eq $transcript_id, @transcript_ids) {
		    push(@transcript_ids,$transcript_id);
			$data_store{'mRNAs'} += 1;
		}
		
		my @gvalues = $feature->get_tag_values('gene_id');
		my $gene_id = shift @gvalues ;
		
		unless (my ($matched) = grep $_ eq $gene_id, @gene_ids) {
		    push(@gene_ids,$gene_id);
			$data_store{'genes'} += 1;
		}
		
	}
}

$gff_file->close();

print "#genes: $data_store{'genes'}\n";
print "#transcripts: $data_store{'mRNAs'}\n";
print "#exons: $data_store{'exons'}\n" ;

# --------------

sub msg {
  my $t = localtime;
  my $line = "[".$t->hms."] @_\n";
  print LOG $line if openhandle(\*LOG);
}

sub runcmd {
  msg("Running:", @_);
  system(@_)==0 or err("Could not run command:", @_);
}

sub err {
  msg(@_);
  exit(2);
}


