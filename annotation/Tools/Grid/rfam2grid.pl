#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Scalar::Util qw(openhandle);
use Time::Piece;
use Time::Seconds;
use FindBin;
use lib ("$FindBin::Bin/PerlLib", "$FindBin::Bin/PerlLibAdaptors");
use File::Basename;
use Bio::SeqIO;
use Cwd;
use Carp;
use Bio::SeqFeature::Generic;
use Bio::Tools::GFF;
no strict qw(subs refs);


my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--fasta filename]
		The name of the genome file to read. 

  Ouput:    
    [--outdir name]
        The name of the output directory. 
		
};

my $grid_computing_module = "BilsGridRunner";
my $rfam_cm_file = "/projects/references/databases/rfam/11.0/models_1_1/E_plus.cm"; # Eukaryote ncRNA models plus LUCA ncRNAs and excluding tRNA (to be annotated by tRNAscan). 
my $gff_formatter = Bio::Tools::GFF->new(-gff_version => 3);

my $outdir = undef;
my $fasta = undef;
my @cmds = ();				# Stores the commands to send to farm
my $quiet;
my @annotations = ();		# Stores Rfama annotations as hashes
my $help;

GetOptions(
    "help" => \$help,
    "fasta=s" => \$fasta,
    "outdir=s" => \$outdir);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

# .. Check that all binaries are available in $PATH

my @tools = ("cmsearch" );	# List of tools to check for!
foreach my $exe (@tools) { check_bin($exe) == 1 or die "Missing executable $exe in PATH"; }

# .. Create output directory 

if (-d $outdir ) {
	die "Output directory $outdir exists. Please remove and try again";
} else {
	msg("Creating output directory $outdir");
	runcmd("mkdir -p $outdir")
}

# .. set up log file

my $logfile = "$outdir/rfam_search.log";
msg("Writing log to: $logfile");
open LOG, '>', $logfile or err("Can't open logfile");

# .. load grid module (courtesy of Brian Haas)

my $perl_lib_repo = "$FindBin::Bin/../PerlLibAdaptors";
msg("-importing module: $grid_computing_module\n");
require "$grid_computing_module.pm" or die "Error, could not import perl module at run-time: $grid_computing_module";

my $grid_computing_method = $grid_computing_module . "::run_on_grid" or die "Failed to initialize GRID module\n";

# .. Read genome fasta file.
my $inseq = Bio::SeqIO->new(-file   => "<$fasta", -format => 'fasta');

# .. and create chunks
msg("Creating chunks for grid\n");

my $seq;

my $seq_counter = 0;

while( $seq = $inseq->next_seq() ) {
	$seq_counter += 1;		
	my $outfile = $outdir . "/seq_" . $seq_counter . ".fasta" ; # We could also use the display_id, but this can cause trouble with special characters
	my $seq_out = Bio::SeqIO->new(-file => ">$outfile" , -format => 'fasta');
	$seq_out->write_seq($seq);
	my $command = "cmsearch --cpu 1 --rfam --cut_tc --tblout " . $outfile . ".rfam "  . $rfam_cm_file . " " . $outfile . " > /dev/null" ;
	push(@cmds,$command);
}

# Submit job chunks to grid

msg("Sending $seq_counter jobs to LSF grid\n");

chomp(@cmds); # Remove empty indices
&$grid_computing_method(@cmds);

# ..Postprocessing here, merging of output and printing gff

msg("Merging output and writing GFF file");

my @files = <$outdir/*.rfam>; 

foreach my $file (@files) {

	open (my $IN, '<', $file) or die "FATAL: Can't open file: $file for reading.\n$!\n";

	while (<$IN>) {
		chomp; 
		my $line = $_; 
		next if ($line =~ /^#.*$/); # Skipping comment lines
		
		my $annotation = parse_line($line);
		push(@annotations,$annotation);
		
	}		
}

my $outfile = "rfam.gff";
open (my $OUT, '>', $outfile) or die "FATAL: Can't open file: $outfile for reading.\n$!\n";

foreach my $feature (@annotations) {
	$feature->gff_format($gff_formatter);
	print $OUT $feature->gff_string, "\n";
}

close($OUT);

# --------------------

sub parse_line {
	# chomp;
	my $line = shift ;
	my ($tn,$tacc,$qn,$qacc,$mdl,$mdlf,$mdlt,$seqf,$seqt,$strand,$trunc,$pass,$gc,$bias,$score,$evalue,$inc,$desc) = split(/\s+/,$line);
	
	my %tags = ( 'rfam-id' => $qn,
			 'rfam-acc' => ($qacc || 'unknown'),
			 'model_start' => $mdlf,
			 'model_end' => $mdlf,
			 'gc-content' => $gc,
			 'ID' => $qacc . "_" .  $tn . "_" . $seqf,
			 'Name' => $qacc . "_" .  $tn . "_" . $seqf,
	    );

		my($from,$to) = sort($seqf,$seqt); # cmsearch reports coordinates in orientation of annotation, not chromosome. Need to sort from low to high for gff
		
	    if( $evalue =~ /[0-9]/ ) {
		$tags{'evalue'} = $evalue;
	    }

	    my $f = Bio::SeqFeature::Generic->new( -seq_id => $tn,
						   -start => $from,
						   -end => $to,
						   -strand => $strand,
						   -frame => 0,
						   -primary_tag => 'ncRNA',	# may argue over whether this is an exon feature, but anything else will be ignored by Maker
						   -source_tag => 'Rfam',
						   -score => $score,
						   -tag => \%tags,
		);
	
		return $f;
}

# --------------------

sub msg {
  my $t = localtime;
  my $line = "[".$t->hms."] @_\n";
  print LOG $line if openhandle(\*LOG);
  print STDERR $line unless $quiet;
}

# --------------------

sub runcmd {
  msg("Running:", @_);
  system(@_)==0 or err("Could not run command:", @_);
}

# --------------------

sub check_bin {
	length(`which @_`) > 0 ? return 1 : return 0;
}

#----------------------------------------------------------------------

sub err {
  $quiet=0;
  msg(@_);
  exit(2);
}
