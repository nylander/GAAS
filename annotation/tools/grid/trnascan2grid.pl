#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Scalar::Util qw(openhandle);
use Time::Piece;
use Time::Seconds;
use FindBin;
use File::Basename;
use Bio::SeqIO;
use Cwd;
use Carp;
use Bio::SeqFeature::Generic;
use Bio::Tools::GFF;
no strict qw(subs refs);
use GAAS::Grid::Bsub;
use GAAS::Grid::Sbatch;
use GAAS::GAAS;

my $header = get_gaas_header();
my $gff_formatter = Bio::Tools::GFF->new(-gff_version => 3);
my $outdir = undef;
my $fasta = undef;
my @cmds = ();				# Stores the commands to send to farm
my $quiet;
my $grid="Slurm";
my @annotations = ();		# Stores Rfama annotations as hashes
my $help;
my $queue=undef;

if ( !GetOptions(
    "h|help!" => \$help,
    "fasta|f=s" => \$fasta,
    "grid=s"  => \$grid,
		"queue=s"  => \$queue,
		"quiet|q!"  => \$quiet,
    "outdir|o=s" => \$outdir))

{
    pod2usage( { -message => "Failed to parse command line",
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($help) {
        pod2usage( { -verbose => 99,
                 -exitval => 0,
                 -message => "$header \n" } );
}

if ( ! (defined($fasta) and defined($outdir) ) ){
    pod2usage( {
           -message => "$header\nAt least 2 parameter are mandatory:\nInput fasta file and output directory \n\n",
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

my @tools = ( "tRNAscan-SE" );	# List of tools to check for!
foreach my $exe (@tools) { check_bin($exe) == 1 or die "Missing executable $exe in PATH"; }

# .. Create output directory

if (-d $outdir ) {
	die "Output directory $outdir exists. Please remove and try again";
} else {
	msg("Creating output directory $outdir");
	runcmd("mkdir -p $outdir")
}

# .. set up log file

my $logfile = "$outdir/trna_search.log";
msg("Writing log to: $logfile");
open LOG, '>', $logfile or err("Can't open logfile");

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
	my $command = "tRNAscan-SE -o " . $outfile . ".trna -q $outfile > /dev/null" ;
	push(@cmds,$command);
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

# ..Postprocessing here, merging of output and printing gff

msg("Merging output and writing GFF file");

my @files = <$outdir/*.trna>;


foreach my $file (@files) {

	open (my $IN, '<', $file) or die "FATAL: Can't open file: $file for reading.\n$!\n";

	while (<$IN>) {
		chomp;
		my $line = $_;
		next if ($line =~ /^Sequence.*/ or $line =~ /^Name.*/ or $line =~ /^---.*/); # Skipping comment lines

		my $annotation = parse_line($line);
		push(@annotations,$annotation);

	}
}

my $outfile = $outdir . "/trnascan_annotations.gff";
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
	my ($seq,$trna_num,$tstart,$tend,$ttype,$tanti,$istart,$iend,$score) = split(/\s+/,$line);

	my %tags = ( 'tRNA-type' => $ttype,
			 'anti-codon' => ($tanti || 'unknown'),
			 'ID' => 'tRNA_' . $ttype . "_" . $seq . "_" . $tstart,
			 'Name' => 'tRNA_' . $ttype . "_" . $seq . "_" . $tstart,
	    );

	my($from,$to) = sort($tstart,$tend);

	my $strand = ( $tstart < $tend ? '1' : '-1' );

	my $f = Bio::SeqFeature::Generic->new( -seq_id => $seq,
						   -start => $from,
						   -end => $to,
						   -strand => $strand,
						   -frame => 0,
						   -primary_tag => 'tRNA',
						   -source_tag => 'tRNAscan',
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

__END__

=head1 NAME

gaas_trnascan2grid.pl

=head1 DESCRIPTION

Chunk input data to run multiple trnascan jobs in parallel using Grid

=head1 SYNOPSIS

    gaas_trnascan2grid.pl -f fasta_file -o outdir
    gaas_trnascan2grid.pl --help

=head1 OPTIONS

=over 8

=item B<--fasta> or B<-f>

The name of the protein fasta file to read.

=item B<--queue>

If you want to define a particular queue to run the jobs

=item B<--grid>

Define which grid to use, Slurm, Lsf or None. Default = Slurm.

=item B<--queue>

If you want to define a particular queue to run the jobs

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
