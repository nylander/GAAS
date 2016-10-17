#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Scalar::Util qw(openhandle);
use Time::Piece;
use Time::Seconds;
use FindBin;
use lib ( "$FindBin::Bin/PerlLib", "$FindBin::Bin/PerlLibAdaptors" );
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
    [--fasta filename]
        The name of the protein file to read.
    [--db name]
        The name of the database use to blast.
    [--chunk_size]
        The number of sequence by job.  If not provided, default size
        will be 500.
    [--nb_seq]
        The number of proteins contained in the db. Useful to cheat on
        the database size. (OrthoMCL aggregation as example). If not
        provided, the current database size is used.
    [--eval]
        The evalue of the sequences kept in the result
Ouput:
    [--outdir name]
        The name of the output directory.

};

my $scipio_outfile = "scipio.merged.gff";
my $outdir         = undef;
my $db             = undef;
my $fasta          = undef;
my $chunk_size     = 500;
my $eval           = 1e-5;
my $nb_seq         = undef;              # Partition size of fasta input
my @chunks = ();    # Holds chunks, partitioning the fasta input (so we
                    # don't send 50.000 jobs to the farm...
my @cmds   = ();    # Stores the commands to send to farm
my $quiet;
my $help;

GetOptions( "help"         => \$help,
            "fasta=s"      => \$fasta,
            "db=s"         => \$db,
            "chunk_size=i" => \$chunk_size,
            "nb_seq=i"     => \$nb_seq,
            "eval"         => \$eval,
            "outdir=s"     => \$outdir );

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}
#check all parameters
if ( !( defined($db) ) ) {
    pod2usage( { -message => "\nA database must be provided (--db)\n",
                 -verbose => 0,
                 -exitval => 2 } );
}
if ( !( defined($fasta) ) ) {
    pod2usage( { -message => "\nA fasta file must be provided (--fasta)\n",
                 -verbose => 0,
                 -exitval => 2 } );
}
if ( !( defined($outdir) ) ) {
    pod2usage( { -message => "\nA database must be provided (--outdir)\n",
                 -verbose => 0,
                 -exitval => 2 } );
}

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

# .. load grid module (courtesy of Brian Haas)

my $grid_computing_module = "BilsGridRunner";
my $perl_lib_repo         = "$FindBin::Bin/../PerlLibAdaptors";
msg("-importing module: $grid_computing_module\n");
require "$grid_computing_module.pm" or
  die
"Error, could not import perl module at run-time: $grid_computing_module";

my $grid_computing_method = $grid_computing_module . "::run_on_grid" or
  die "Failed to initialize GRID module\n";

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
chomp(@cmds);    # Remove empty indices
&$grid_computing_method(@cmds);

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
