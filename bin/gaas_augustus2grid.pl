#!/usr/bin/env perl

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
no strict qw(subs refs);


my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
	[--hints filename ]
		Intron hints for Augustus
	[--species name ]
		Species name compatible with Augustus
	[--genome filename]
		The name of the genome file to align to. 
  Ouput:    
    [--outdir name]
        The name of the output directory. 
		
};

my $outdir = undef;
my $genome = undef;
my $species = undef;
my $hints = undef;
my @cmds = ();				# Stores the commands to send to farm
my $quiet;
my $help;

GetOptions(
    "help" => \$help,
    "hints=s" => \$hints,
	"species=s" => \$species,
	"genome=s" => \$genome,
    "outdir=s" => \$outdir);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

# .. Check that all binaries are available in $PATH

my @tools = ( "augustus" );	# List of tools to check for!
foreach my $exe (@tools) { check_bin($exe) == 1 or die "Missing executable $exe in PATH"; }

my $augustus_config_pathj = $ENV{'AUGUSTUS_CONFIG_PATH'} or die "AUGUSTUS_CONFIG_PATH is not set, aborting." ;

# .. Create output directory 

if (-d $outdir ) {
	die "Output directory $outdir exists. Please remove and try again";
} else {
	msg("Creating output directory $outdir");
	runcmd("mkdir -p $outdir")
}

# .. set up log file

my $logfile = "$outdir/augustus.log";
msg("Writing log to: $logfile");
open LOG, '>', $logfile or err("Can't open logfile");

# .. load grid module (courtesy of Brian Haas)

my $grid_computing_module = "BilsGridRunner";
my $perl_lib_repo = "$FindBin::Bin/../PerlLibAdaptors";
msg("-importing module: $grid_computing_module\n");
require "$grid_computing_module.pm" or die "Error, could not import perl module at run-time: $grid_computing_module";

my $grid_computing_method = $grid_computing_module . "::run_on_grid" or die "Failed to initialize GRID module\n";

# .. Read genome fasta file.
my $inseq = Bio::SeqIO->new(-file   => "<$genome", -format => 'fasta');

# .. and create chunks
msg("Creating chunks for GRID\n");

my $counter = 10000;
my $seq;

while( $seq = $inseq->next_seq() ) {
	$counter += 1;
    my $outfile = $outdir . "/seq_" . $counter . ".fa";
	my $cmd = "augustus --species=$species --hintsfile=$hints --alternatives-from-evidence=true --gff3=on --extrinsicCfgFile=/references/software/augustus/config/extrinsic/extrinsic.E.cfg --uniqueGeneId=true $outfile > $outfile.augustus" ;
	push(@cmds,$cmd);
	my $seq_out = Bio::SeqIO->new(-file => ">$outfile", -format => 'fasta');
	$seq_out->write_seq($seq);
}

# Submit job chunks to grid

chomp(@cmds); # Remove empty indices
&$grid_computing_method(@cmds);


# ..Postprocessing here, like merging of output files

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
