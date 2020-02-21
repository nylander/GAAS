#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Bio::DB::EUtilities;
use Bio::SeqIO;
use Scalar::Util qw(openhandle);
use Time::Piece;
use Time::Seconds;
use Pod::Usage;


my $header = qq{
########################################################
# NBIS 2015 - Sweden                                   #
# Marc Hoeppner / Jacques Dainat                        #
# Please cite NBIS (www.nbis.se) when using this tool. #
########################################################
};

my $outfile = undef;
my $format = "fasta";
my $quiet;
my $organisms = undef;
my $dbs = undef;
my $outdir = "tmp";
my @dbs = ( "nucest" , "protein" );

my $list;
my $help;

if ( !GetOptions(
    "help|h" => \$help,
	"l|list" => \$list,
	"outdir=s" => \$outdir,
	"f|format=s" => \$format,
	"o|organisms=s" => \$organisms,
	"db|dbs=s"	=> \$dbs,
    "o|output|outfile=s" => \$outfile))
{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($help) {
    pod2usage( { -verbose => 2,
                 -exitval => 2,
                 -message => "$header\n" } );
}


if ($list) {
	my $db_factory = Bio::DB::EUtilities->new(-eutil => 'einfo', -email => 'mymail@foo.bar',);
	my $db_list = join("\n\t",$db_factory->get_available_databases);
	print "\n", $db_list , "\n";
	exit(0);
}

if (! defined($organisms) ){
    pod2usage( {
           -message => "$header\nAt least 1 parameter is mandatory:\nInput organisms (-o species1:species2:species3).\n\n".
           "Ouptut is optional. Look at the help documentation to know more.\n",
           -verbose => 0,
           -exitval => 1 } );
}

# .. Create output directory

runcmd("mkdir -p $outdir");

# .. set up log file

my $logfile = "$outdir/reference_sequences.log";
msg("Writing log to: $logfile");
open LOG, '>', $logfile or err("Can't open logfile");

if ($dbs) {
	@dbs = split(":",$dbs);
}

# Iterate over all organisms (can be one..)
foreach my $organism (split(":", $organisms)) {

	my $query_term = join(" ",split("_",$organism)) . "[ORGN]";

	foreach my $db (@dbs) {

		my $factory = Bio::DB::EUtilities->new(-eutil      => 'esearch',
		                                       -email      => 'me@foo.com',
		                                       -db         => $db,
											   -retmax 	   => [10],
		                                       -term       => $query_term,
		                                       -usehistory => 'y');

		my $count = $factory->get_count;

		msg("Found " . $factory->get_count . " hits for " . $organism . " in database '" . $db . "'\n");

		next if ($count == 0); # Skip if nothing was found

		my $hist  = $factory->next_History || die 'No history data returned';
		print "History returned\n";

		# note db carries over from above
		$factory->set_parameters(-eutil   => 'efetch',
		                         -rettype => $format,
		                         -history => $hist);

		my $retry = 0;
		my ($retmax, $retstart) = (500,0);

		open (my $out, '>', $organism . "." . $db . ".fa") || die "Can't open file:$!";

		RETRIEVE_SEQS:

		while ($retstart < $count) {
		    $factory->set_parameters(-retmax   => $retmax,
			                         -rettype => $format,
									 -retstart => $retstart);
		    eval{
		        $factory->get_Response(-cb => sub {my ($data) = @_; print $out $data} );
		    };
		    if ($@) {
		        die "Server error: $@.  Try again later" if $retry == 5;
		        print STDERR "Server error, redo #$retry\n";
		        $retry++ && redo RETRIEVE_SEQS;
		    }
		    print "Retrieved $retstart";
		    $retstart += $retmax;
		}
		close $out;
		print "\n";


	} # end databases

} # end organism

#######################################################################################################################
        ####################
         #     METHODS    #
          ################
           ##############
            ############
             ##########
              ########
               ######
                ####
                 ##
sub msg {
  my $t = localtime;
  my $line = "[".$t->hms."] @_\n";
  print LOG $line if openhandle(\*LOG);
  print STDERR $line unless $quiet;
}

sub runcmd {
  msg("Running:", @_);
  system(@_)==0 or err("Could not run command:", @_);
}

=head1 NAME

ncbi_get_reference_data.pl -
The script allow to recovered information from NCBI databases.
The result is written to the specified output file, or to STDOUT.

=head1 SYNOPSIS

    ./ncbi_get_reference_data.pl -o species1:species2:species3 [ -o outfile ]
    ./ncbi_get_reference_data.pl --help

=head1 OPTIONS

=over 8

=item B<-l> or B<--list>

List of all available databases

=item B<-o> or B<--organisms>

The names of the species to query data from. Species name format: Genus_species (e.g. Gallus_gallus). When querying several organisms please follow this nomenclature: species1:species2:species3

=item B<--db> or B<--dbs>

The names of the NCBI databases to query for data. Default: nucest, protein (see --list for options). When querying several databases please follow this nomenclature: db1:db2:db3

=item B<-f> or B<--format>

The file format to produce. Not all databases can write all formats! Default: fasta

=item B<-o>, B<--output> or B<--outfile>

The name of the output file. By default the output is the standard output

=item B<-h> or B<--help>

Display this helpful text.

=back

=cut
