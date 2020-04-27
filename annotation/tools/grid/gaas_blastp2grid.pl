#!/usr/bin/env perl

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Scalar::Util qw(openhandle);
use Time::Piece;
use Time::Seconds;
use File::Basename;
use Bio::SeqIO;
use Cwd;
use Carp;
no strict qw(subs refs);
use GAAS::Grid::Bsub;
use GAAS::Grid::Sbatch;
use GAAS::GAAS;

my $header = get_gaas_header();
my $outdir         = "blastp_output";
my $db             = undef;
my $fasta          = undef;
my $chunk_size     = 500;
my $eval           = 1e-5;
my $nb_seq         = undef;              # Partition size of fasta input
my @chunks = ();    # Holds chunks, partitioning the fasta input (so we
                    # don't send 50.000 jobs to the farm...
my @cmds   = ();    # Stores the commands to send to farm
my $quiet;
my $grid="Slurm";
my $queue=undef;
my $help;

if ( !	GetOptions( "h|help!"      => \$help,
				            "fasta|f=s"      => \$fasta,
				            "db=s"         => \$db,
				            "chunk_size=i" => \$chunk_size,
				            "nb_seq=i"     => \$nb_seq,
				            "eval"         => \$eval,
										"grid=s"  		 => \$grid,
										"quiet|q!" 		 => \$quiet,
										"queue=s"  		 => \$queue,
				            "outdir=s"     => \$outdir ) )

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

if ( ! defined($db) or  ! defined($fasta) ){
    pod2usage( {
           -message => "$header\nAt least 2 parameters are mandatory:\nInput fasta file (--fasta) and a database (--db)\n\n",
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

my @tools = ("blastp");
foreach my $exe (@tools) {
    check_bin($exe) == 1 or die "Missing executable $exe in PATH";
}

# .. Create output directory

if ( -d $outdir ) {
    msg( "Be careful, we are using an existinf Output directory $outdir. " .
         "If you do not want that, you have to stop the job" );
}
else {
    msg("Creating output directory $outdir");
    runcmd("mkdir -p $outdir");
}

# .. set up log file

my $logfile = "$outdir/blastp2grid.log";
msg("Writing log to: $logfile");
open LOG, '>', $logfile or err ("Can't open logfile");

# .. Read protein fasta file.
my $inseq = Bio::SeqIO->new( -file => "<$fasta", -format => 'fasta' );

# .. and create chunks
msg("Creating chunks for GRID\n");

my @seqarray      = ();
my $counter       = 0;
my $chunk_counter = 1;

my $seq;

while ( $seq = $inseq->next_seq() ) {
    $counter += 1;
    push( @seqarray, $seq );

    if ( $counter == $chunk_size ) {
        my $outfile = $outdir . "/chunk_" . $chunk_counter . ".fa";
        write_chunk( $outfile, @seqarray );
        @seqarray = ();
        $chunk_counter += 1;
        $counter = 0;
    }
}
my $outfile =
  $outdir . "/chunk_" .
  $chunk_counter . ".fa";   # Clunky, the last chunk is <= chunk_size...
write_chunk( $outfile, @seqarray );

# Push all jobs into the command list
if ( !defined($nb_seq) ) {
    for ( my $i = 1; $i <= $chunk_counter; $i++ ) {
        my $cmd =
          "blastp -evalue $eval -num_alignments 100000 " .
          "-seg yes -outfmt 6 -db $db -query $outdir/chunk_$i.fa " .
          "-out $outdir/chunk_$i.tab";
        push( @cmds, $cmd );
    }
}
else {
    for ( my $i = 1; $i <= $chunk_counter; $i++ ) {
        my $cmd =
          "blastp -dbsize $nb_seq -evalue $eval " .
          "-num_alignments 100000 -seg yes -outfmt 6 -db $db " .
          "-query $outdir/chunk_$i.fa -out $outdir/chunk_$i.tab";
        push( @cmds, $cmd );
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

# Merging the outputs
msg("Merging outputs from chunks");

my @files = <$outdir/*.tab>;

foreach my $file (@files) {
    system("cat $file >> $outdir/blastp.merged");
}

msg("Finished BLASTp grid run.");

# --------------------

sub write_chunk
{
    my $outfile = shift;
    my @seqs    = @_;
    my $seq_out =
      Bio::SeqIO->new( -file => ">$outfile", -format => 'fasta' );
    foreach my $seq (@seqs) { $seq_out->write_seq($seq) }
}

# --------------------

sub msg
{
    my $t    = localtime;
    my $line = "[" . $t->hms . "] @_\n";
    print LOG $line if openhandle( \*LOG );
    print STDERR $line unless $quiet;
}

# --------------------

sub runcmd
{
    msg( "Running:", @_ );
    system(@_) == 0 or err ( "Could not run command:", @_ );
}

# --------------------

sub check_bin
{
    length(`which @_`) > 0 ? return 1 : return 0;
}

#----------------------------------------------------------------------

sub err
{
    $quiet = 0;
    msg(@_);
    exit(2);

}

__END__

=head1 NAME

gaas_blastp2grid.pl -

=head1 DESCRIPTION

Chunk input data to run multiple blastp jobs in parallel to grid

=head1 SYNOPSIS

    gaas_blastp2grid.pl -f fasta_file --db db_name
    gaas_blastp2grid.pl --help

=head1 OPTIONS

=over 8

=item B<--fasta> or B<-f>

The name of the protein fasta file to use as query.

=item B<--db>

The name of the database use to blast

=item B<--nb_seq>

The number of proteins contained in the db. Useful to cheat on
the database size. (OrthoMCL aggregation as example). If not
provided, the current database size is used.

=item B<--eval>

The maximu evalue of the sequences kept in the result

=item B<--chunk_size>

The number of sequence by job. If not provided, default size
will be 500.

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
