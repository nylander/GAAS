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
use Cwd;
use Carp;
no strict qw(subs refs);


my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]

  Input:
    [--query filename]
	The name of the query genome file.
    [--target filename]
	The name of the target genome file. 

  Ouput:    
    [--outdir name]
        The name of the output directory. 
		
};

my $grid_computing_module = "BilsGridRunner";

my $outdir = undef;
my $query = undef;
my $target = undef;
my @cmds = ();				# Stores the commands to send to farm
my @query_seqs = ();			# List of sequences from the query genome
my @target_seqs = ();			# List of sequences from the target genome
my @lav_files = ();
my $job_limit = 500;			# Maximum number of jobs to allow before aborting
my $quiet;
my $help;

GetOptions(
    "help" => \$help,
    "query=s" => \$query,
    "target=s" => \$target,
    "outdir=s" => \$outdir);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

# .. Check that all binaries are available in $PATH

my @tools = ("faToNib" , "lastz" , "fastaexplode" );	# List of tools to check for!
foreach my $exe (@tools) { check_bin($exe) == 1 or die "Missing executable $exe in PATH"; }

my $working_dir = getcwd;

my $query_base = basename($query) ;
my $target_base = basename($target) ;
 
my $query_dir = $outdir . "/" . $query_base ;
my $target_dir = $outdir . "/" . $target_base ;


# .. Create output directory 

if (-d $outdir ) {
	die "Output directory $outdir exists. Please remove and try again";
} else {
	msg("Creating output directory $outdir");
	runcmd("mkdir -p $outdir");
	runcmd("mkdir -p $outdir/chain");
	runcmd("mkdir -p $outdir/lav");
	runcmd("mkdir -p $outdir/psl");
	runcmd("mkdir -p $query_dir");
	runcmd("mkdir -p $target_dir");
}

# .. set up log file

my $logfile = "$outdir/lastz.log";
msg("Writing log to: $logfile");
open LOG, '>', $logfile or err("Can't open logfile");

# .. load grid module (courtesy of Brian Haas)

my $perl_lib_repo = "$FindBin::Bin/../PerlLibAdaptors";
msg("-importing module: $grid_computing_module\n");
require "$grid_computing_module.pm" or die "Error, could not import perl module at run-time: $grid_computing_module";

my $grid_computing_method = $grid_computing_module . "::run_on_grid" or die "Failed to initialize GRID module\n";

msg("Generating size indices for genomes");

my $query_size = $outdir . "/" . basename($query) . ".sizes" ;
my $target_size = $outdir . "/" . basename($target) . ".sizes" ; 

runcmd("faSize $query -detailed > $query_size");
runcmd("faSize $target -detailed > $target_size");

msg("Exploding genome sequences for recursive searches");

runcmd("fastaexplode -f $query -d $query_dir");
runcmd("fastaexplode -f $target -d $target_dir");

msg("Converting and gathering sequences for LastZ recursive alignments...");

msg("Searching $query_dir");

opendir (DIR, $query_dir ."/" ) or die "Could not open directory for reading ($query_dir)" ;

while (my $file = readdir(DIR)) {
	next unless ($file =~ m/fa$/ );
	my $nib_file = $file;
        $nib_file =~ s/fa$/nib/ ;
        system("faToNib $query_dir/$file $query_dir/$nib_file");
        push @query_seqs, $nib_file ;

}

closedir(DIR);

opendir(DIR, $target_dir) or die "Could not open directpry for reading ($target_dir)";

msg("Searching $target_dir");

while (my $file = readdir(DIR)) {
        next unless ($file =~ m/fa$/ );
        my $nib_file = $file;
        $nib_file =~ s/fa$/nib/ ;
	system("faToNib $target_dir/$file $target_dir/$nib_file");
        push @target_seqs, $nib_file ;
}

closedir(DIR);

my $query_jobs = scalar @query_seqs ;
my $target_jobs = scalar @target_seqs ;

die "Way too many jobs - can't submit that to the grid. Consider limiting your input data" if ($target_jobs * $query_jobs > $job_limit);

msg("Building commands for GRID...");

foreach my $query_seq(@query_seqs) {

	foreach my $target_seq(@target_seqs) {
		my $lav_file = $query_seq . "-" . $target_seq . ".lav" ;
		push @lav_files , $lav_file ; 
		push @cmds , "lastz $query_dir/$query_seq $target_dir/$target_seq > $outdir/lav/$lav_file" ;
	}	

}


# Submit job chunks to grid

msg("Sending jobs to LSF grid\n");

chomp(@cmds); # Remove empty indices

&$grid_computing_method(@cmds);

msg("### Converting LAV files to PSL ###");

foreach my $lav_file (@lav_files) {

	my $psl_file = $lav_file ;
	$psl_file =~ s/lav$/psl/ ;  
	runcmd("lavToPsl $outdir/lav/$lav_file $outdir/psl/$psl_file");
}

opendir(DIR, $outdir."/psl") or die "Could not open directory for reading ($outdir)" ;


while (my $file = readdir(DIR)) {

	next unless ($file =~ m/psl$/ );

	my $chainfile = $file;
	$chainfile =~ s/psl$/chain/ ;

	runcmd("axtChain -linearGap=loose -psl $outdir/psl/$file $query_dir $target_dir $outdir/chain/$chainfile");

}

closedir(DIR);

# Merge Chains

msg("### Merging chains ###");

runcmd("chainMergeSort $outdir/chain/*.chain > $outdir/all.chain");

runcmd("chainPreNet $outdir/all.chain $query_size $target_size $outdir/all.pre.chain");

# Netting chains

msg("### Netting the chains ###");

runcmd("chainNet $outdir/all.pre.chain -minSpace=1 $query_size $target_size stdout /dev/null | netSyntenic stdin $outdir/lastz.net");



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
