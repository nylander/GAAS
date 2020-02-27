#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Scalar::Util qw(openhandle);
use Time::Piece;
use POSIX;
use Time::Seconds;
use NBIS::Grid::Bsub;
use NBIS::Grid::Sbatch;
use File::Basename;
use Bio::SeqIO;
use Cwd;
use Carp;
use Bio::SeqFeature::Generic;
no strict qw(subs refs);

my $header = qq{
########################################################
# NBIS - Sweden                                        #
#                                                      #
# Please cite NBIS (www.NBIS.se) when using this tool. #
########################################################
};

my $outdir = transposonPSI_output;
my $fastaFile;
my @cmds = ();				# Stores the commands to send to farm
my $quiet;
my $help;
my $grid="Slurm";
my $queue=undef;
my $chunk=10;

if ( !GetOptions(
    "help" => \$help,
    "fasta|f=s" => \$fastaFile,
    "grid=s"  => \$grid,
    "chunk=i"  => \$chunk,
    "queue=s"  => \$queue,
    "outdir|o=s" => \$outdir))

{
    pod2usage( { -message => "Failed to parse command line",
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($help) {
        pod2usage( { -verbose => 1,
                 -exitval => 0,
                 -message => "$header \n" } );
}

if ( ! defined($fastaFile) ){
    pod2usage( {
           -message => "$header\nAt least 2 parameter are mandatory:\nInput fasta file and output directory \n\n",
           -verbose => 0,
           -exitval => 2 } );
}

if (! -e $fastaFile){
	print "The fasta file ".$fastaFile." does not exist.\n";exit;
}

# set grid
my @grid_choice=('slurm','lsf','none');
$grid=lc($grid);
if (! grep( /^$grid/, @grid_choice ) ) {
  print "$grid is not a value accepted for grid parameter.";exit;
}
$grid= undef if lc($grid) eq 'none';

# .. Check that all binaries are available in $PATH
my @tools = ("transposonPSI.pl" );	# List of tools to check for!
foreach my $exe (@tools) { check_bin($exe) == 1 or die "Missing executable $exe in PATH"; }

# .. Create output directory

if (-d $outdir ) {
	die "Output directory $outdir exists. Please remove and try again";
} else {
	msg("Creating output directory $outdir");
	runcmd("mkdir -p $outdir")
}

# .. set up log file

my $logfile = "$outdir/transposonPSI.log";
msg("Writing log to: $logfile");
open LOG, '>', $logfile or err("Can't open logfile");

# .. Read genome fasta file.
my $inseq = Bio::SeqIO->new(-file   => "<$fastaFile", -format => 'fasta');

# .. and create chunks
msg("Creating chunks\n");

my $seq;

my $seq_counter = 0;
my $nbLine=`grep -c ">" $fastaFile`;
my $chunk_size= ceil($nbLine/$chunk);
my $chunk_counter=0;
my $seq_out;
my $chunk_file ;

while( $seq = $inseq->next_seq() ) {
	if ($chunk_counter == 0){
		$seq_counter += 1;
		$chunk_file = $outdir."/".$fastaFile. "_chunck" . $seq_counter . ".fasta" ; # We could also use the display_id, but this can cause trouble with special characters
		$seq_out = Bio::SeqIO->new(-file => ">$chunk_file" , -format => 'fasta');
	}

	$seq_out->write_seq($seq);
	$chunk_counter++;

	if ($chunk_counter == $chunk_size){
		my $command = "transposonPSI.pl ".$chunk_file." prot > /dev/null" ;
		push(@cmds,$command);
		$chunk_counter = 0;
	}
}

msg("submitting chunks\n");

if( $grid){
  msg("Sending $seq_counter jobs to the grid\n");
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

# ..Postprocessing here, merging of output and printing

msg("Merging outputS");

# Deal with TPSI.allHits
my $outfile = $fastaFile.".all.TPSI.allHits";
my @files_allHits = <*.TPSI.allHits>;
open (my $OUT, '>', $outfile) or die "FATAL: Can't open file: $outfile for reading.\n$!\n";

foreach my $file (@files_allHits) {

	open (my $IN, '<', $file) or die "FATAL: Can't open file: $file for reading.\n$!\n";

	while (<$IN>) {
		my $line = $_;
		next if ($line =~ /^\/\/.*$/); # Skipping comment lines
		print $OUT $line;
	}
}
close($OUT);
my $command = "mv *.TPSI.allHits $outdir";
system ("/bin/bash -c '$command'");

# Deal with TPSI.topHits
$outfile = $fastaFile.".all.TPSI.topHits";
my @files_topHits = <*.TPSI.topHits>;
open ($OUT, '>', $outfile) or die "FATAL: Can't open file: $outfile for reading.\n$!\n";

foreach my $file (@files_topHits) {

        open (my $IN, '<', $file) or die "FATAL: Can't open file: $file for reading.\n$!\n";

        while (<$IN>) {
                my $line = $_;
                next if ($line =~ /^\/\/.*$/); # Skipping comment lines
		print $OUT $line;
        }
}
close($OUT);
$command = "mv *.TPSI.topHits $outdir";
system ("/bin/bash -c '$command'");

# --------------------

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

transposonPSI2grid.pl -

=head1 SYNOPSIS

We transposonPSI against protein fasta sequences

    ./transposonPSI2grid.pl -f fasta_file -o outdir
    ./transposonPSI2grid.pl --help

=head1 OPTIONS

=over 8

=item B<--fasta> or B<-f>

The name of the protein fasta file to read.

=item B<--chunk>

By default 10. We slice the fasta input file in many chunk to distribute more efficiently small tasks to each cpu.

=item B<--queue>

If you want to define a particular queue to run the jobs

=item B<--grid>

Define which grid to use, Slurm, Lsf or None. Default = Slurm.

=item B<--outdir> or B<-o>

The name of the output directory.

=item B<-h> or B<--help>

Display this helpful text.

=back

=cut
