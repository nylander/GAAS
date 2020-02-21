#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
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

my $header = qq{
########################################################
# NBIS - Sweden                                        #
#                                                      #
# Please cite NBIS (www.NBIS.se) when using this tool. #
########################################################
};


my $grid_computing_module = "BilsGridRunner";
my $pfam_hmm_file = "/projects/references/databases/pfam/31.0/Pfam-A.hmm";

my $outdir = undef;
my $fasta = undef;
my @cmds = ();				# Stores the commands to send to farm
my $quiet;
my $help;
my $nogrid=undef;
my $chunk_size     = 500;
my @chunks = ();    # Holds chunks, partitioning the fasta input (so we
                    # don't send 50.000 jobs to the farm...

if ( !GetOptions(
    "help" => \$help,
    "fasta|f=s" => \$fasta,
    "hmm=s"  => \$pfam_hmm_file,
    "chunk_size=i" => \$chunk_size,
    "nogrid!"  => \$nogrid,
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

if ( ! (defined($fasta) and defined($outdir) ) ){
    pod2usage( {
           -message => "$header\nAt least 2 parameter are mandatory:\nInput fasta file and output directory \n\n",
           -verbose => 0,
           -exitval => 2 } );
}

if (! -e $pfam_hmm_file){
	print "The cm file ".$pfam_hmm_file." does not exist. Please define it using the cm option.\n";exit;
}
if (! -e $fasta){
	print "The fasta file ".$fasta." does not exist.\n";exit;
}

# .. Check that all binaries are available in $PATH
my @tools = ("hmmscan" );	# List of tools to check for!
foreach my $exe (@tools) { check_bin($exe) == 1 or die "Missing executable $exe in PATH"; }

# .. Create output directory 

if (-d $outdir ) {
	die "Output directory $outdir exists. Please remove and try again";
} else {
	msg("Creating output directory $outdir");
	runcmd("mkdir -p $outdir")
}

# .. set up log file

my $logfile = "$outdir/pfam_search.log";
msg("Writing log to: $logfile");
open LOG, '>', $logfile or err("Can't open logfile");

# .. load grid module (courtesy of Brian Haas)
my $grid_computing_method;
if(! $nogrid){

	my $perl_lib_repo = "$FindBin::Bin/../PerlLibAdaptors";
	msg("-importing module: $grid_computing_module\n");
	require "$grid_computing_module.pm" or die "Error, could not import perl module at run-time: $grid_computing_module";

	$grid_computing_method = $grid_computing_module . "::run_on_grid" or die "Failed to initialize GRID module\n";
}

# .. Read genome fasta file.
my $inseq = Bio::SeqIO->new(-file   => "<$fasta", -format => 'fasta');

# .. and create chunks
msg("Creating chunks\n");

my @seqarray      = ();
my $counter       = 0;
my $chunk_counter = 1;
my $seq;

while( $seq = $inseq->next_seq() ) {
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
for ( my $i = 1; $i <= $chunk_counter; $i++ ) {
     my $infile = $outdir . "/chunk_" . $i . ".fa"; 
     my $outfile = $outdir . "/chunk_" . $i . ".pfam";
     my $cmd = "hmmscan --cpu 1 --domtblout " . $outfile . " "  . $pfam_hmm_file . " " . $infile . " > /dev/null" ;
     push( @cmds, $cmd );
 }

msg("submitting chunks\n");

if( ! $nogrid){
	# Submit job chunks to grid
	msg("Sending $chunk_counter jobs to LSF grid\n");
	chomp(@cmds); # Remove empty indices
	&$grid_computing_method(@cmds);
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

# ..Postprocessing here, merging of output

msg("Merging output and writing GFF file");

my @files = <$outdir/*.pfam>; 

my $outfile = "pfam.merged";
open (my $OUT, '>', $outdir."/".$outfile) or die "FATAL: Can't open file: $outfile for reading.\n$!\n";

foreach my $file (@files) {

        open (my $IN, '<', $file) or die "FATAL: Can't open file: $file for reading.\n$!\n";

        while (<$IN>) {
                chomp; 
                my $line = $_; 
                next if ($line =~ /^#.*$/); # Skipping comment lines
                
		            print $OUT $line;                
        }               
}

close ($OUT);

msg("Finished pfam grid run.");

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

pfam2grid.pl -
We run hmmscan searches against a pfam.hmm 

=head1 SYNOPSIS

    ./pfam2grid.pl -f genome.fasta -o outdir
    ./pfam2grid.pl --help

=head1 OPTIONS

=over 8

=item B<--fasta> or B<-f>

The name of the fasta file to read. 

=item B<--hmm> 

File containing the pfam hmm models 

=item B<--nogrid> 

Do not use the script in grid version.

=item B<--outdir> or B<-o>

The name of the output directory. 

=item B<-h> or B<--help>

Display this helpful text.

=back

=cut
