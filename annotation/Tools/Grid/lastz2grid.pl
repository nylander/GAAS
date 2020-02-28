#!/usr/bin/env perl

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Scalar::Util qw(openhandle);
use Time::Piece;
use Time::Seconds;
use File::Basename;
use Cwd;
use Carp;
no strict qw(subs refs);
use GAAS::Grid::Bsub;
use GAAS::Grid::Sbatch;
use GAAS::GAAS;

my $header = get_gaas_header();
my $outdir = undef;
my $query = undef;
my $target = undef;
my @cmds = ();				# Stores the commands to send to farm
my @query_seqs = ();			# List of sequences from the query genome
my @target_seqs = ();			# List of sequences from the target genome
my @lav_files = ();
my $job_limit = 500;			# Maximum number of jobs to allow before aborting
my $quiet;
my $grid="Slurm";
my $queue=undef;
my $help;

if ( ! GetOptions(
						"query=s" 			=> \$query,
						"target=s" 			=> \$target,
						"grid=s"  			=> \$grid,
						"quiet|q!" 			=> \$quiet,
						"queue=s"  			=> \$queue,
				    "outdir=s" 			=> \$outdir,
						"help|h!"				=> \$help ) )

{
    pod2usage( { -message => "Failed to parse command line",
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($help) {
        pod2usage( { -verbose => 99,
		                 -exitval => 0,
		                 -message => "$header\n" } );
}

if ( ! defined($query) or ! defined($target)){
    pod2usage( {
           -message => "$header\nAt least 2 parameters are mandatory:\n a query genome file (--query) and a target genome file (--target)\n\n",
           -verbose => 0,
           -exitval => 2 } );
}

# set grid option properly
my @grid_choice=('slurm','lsf','none');
$grid=lc($grid);
if (! grep( /^$grid/, @grid_choice ) ) {
  print "$grid is not a value accepted for grid parameter.";exit;
}
$grid= undef if lc($grid) eq 'none';

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
msg("submitting chunks\n");

if( $grid){
	msg("Sending $#cmds jobs to the grid\n");
  chomp(@cmds); # Remove empty indices
  # Submit job chunks to grid
  my $grid_runner;
  if ( $grid eq 'lsf'){
    $grid_runner = Bsub->new( cmds_list => \@cmds);
  }
  elsif( $grid eq 'slurm'){
    $grid_runner = Sbatch->new( cmds_list => \@cmds);
  }
	if($queue){$grid_runner->queue($queue)}
  $grid_runner->run();
}
else{
 	foreach my $command (@cmds){

 		system($command);

 		if ($? == -1) {
    			 print "failed to execute: $!\n";
		}
 		elsif ($? & 127) {
 		    printf "child died with signal %d, %s coredump\n",
 		    ($? & 127),  ($? & 128) ? 'with' : 'without';
 		}
 		else {
 		    printf "child exited with value %d\n", $? >> 8;
 		}
 	}
}

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

__END__

=head1 NAME

gaas_blat2grid.pl

=head1 DESCRIPTION

Chunk input data to run multiple blat jobs in parallel to grid

=head1 SYNOPSIS

    gaas_blat2grid.pl -f fasta_file --db db_name
    gaas_blat2grid.pl --help

=head1 OPTIONS

=over 8

=item B<--query>

The name of the query genome file.

=item B<--target>

The name of the target genome file.

=item B<--queue>

If you want to define a particular queue to run the jobs

=item B<--grid>

Define which grid to use, Slurm, Lsf or None. Default = Slurm.

=item B<--quiet> or B<-q>

Quiet mode

=item B<--outdir> or B<-o>

The name of the output directory.

=item B<-h> or B<--help>

Display this helpful text.

=back

=head1 FEEDBACK

=head2 Did you find a bug?

Do not hesitate to report bugs to help us keep track of the bugs and their
resolution. Please use the GitHub issue tracking system available at this
address:

            https://github.com/NBISweden/GAAS/issues

 Ensure that the bug was not already reported by searching under Issues.
 If you're unable to find an (open) issue addressing the problem, open a new one.
 Try as much as possible to include in the issue when relevant:
 - a clear description,
 - as much relevant information as possible,
 - the command used,
 - a data sample,
 - an explanation of the expected behaviour that is not occurring.

=head2 Do you want to contribute?

You are very welcome, visit this address for the Contributing guidelines:
https://github.com/NBISweden/GAAS/blob/master/CONTRIBUTING.md

=cut

AUTHOR - Jacques Dainat
