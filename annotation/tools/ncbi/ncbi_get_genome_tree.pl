#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Scalar::Util qw(openhandle);
use Time::Piece;
use Time::Seconds;
use Pod::Usage;
use XML::LibXML;
use Data::Dumper;
use LWP::Simple;
use LWP::UserAgent;
use HTTP::Request::Common;
use Try::Tiny;
use Bio::DB::Taxonomy;
use Bio::TreeIO;
use Bio::Tree::Tree;
use Bio::Tree::TreeFunctionsI;
use Bio::DB::EUtilities;
use Bio::SeqIO;
use GAAS::GAAS;

my $header = get_gaas_header();
my $outfile = undef;
my $quiet = undef;
my $message="";
my $taxid=undef;
my $verbose=undef;
my $list;
my $help;

if ( !GetOptions(
    "help|h"    => \$help,
  	"t|taxid=i" => \$taxid,
  	"q"         => \$quiet,
    "v"         => \$verbose,
#	"outdir=s" => \$outdir,
    "o|output|outfile=s" => \$outfile))
{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($help) {
    pod2usage( { -verbose => 99,
                 -exitval => 0,
                 -message => "$header\n" } );
}


# .. Create output directory

#runcmd("mkdir -p $outdir");

# .. set up log file

#my $logfile = "$outdir/reference_sequences.log";
#msg("Writing log to: $logfile");
#open LOG, '>', $logfile or err("Can't open logfile");


##
my %KINGDOM = (1 => 'Eukaryota', 2 => 'Bacteria', 3 => 'Archaea',4 => 'Viroids', 5 => 'Viruses');

my %GROUP = ('Eukaryota' => {1 => 'All', 2 => 'Animals', 3 => 'Fungi', 4 => 'Other', 5 => 'Plants', 6 => 'Protists'},
		'Bacteria' => {1 => 'All', 2 => 'Acidobacteria', 3 => 'Aquificae', 4 => 'Caldiserica', 5 => 'Chrysiogenetes', 6 => 'Deferribacteres', 7 => 'Dictyoglomi', 8 => 'Elusimicrobia', 9 => 'FCB group', 10 => 'Fusobacteira', 11 => 'Nitrospinae/Tectomicrobia group',
		12 => 'Nitrospirae', 13 => 'PVC group', 14 => 'Proteobacteria', 15 => 'Rhodothermaeota', 16 => 'Spirochaetes', 17 => 'Synergistetes', 18 => 'Terrabacteria group', 19 => 'Thermodesulfobacteria', 20 => 'Thermotogae', 21 => 'environmental samples', 22 => 'unclassified Bacteria'},
		'Archaea' => {1 => 'All'},
		'Viroids' => {1 => 'All'},
		'Viruses' => {1 => 'All'});

my %SUBGROUP = ('Animals' => {1 => 'All', 2 => 'Amphibians', 3 => 'Birds', 4 => 'Fishes', 5 => 'Flatworms', 6 => 'Insects', 7 => 'Mammals', 8 => 'Other Animals', 9 => 'Reptiles', 10 => 'Roundworms'},
		'Fungi' => {1 => 'All', 2 => 'Ascomycetes', 3 => 'Basidiomycetes', 4 => 'Other Fungi'},
		'Other' => {1 => 'All'},
		'Plants' => {1 => 'All', 2 => 'Green Algae', 3 => 'Land Plants', 4 => 'Other Plants'},
		'Protists' => {1 => 'All', 2 => 'Apicomplexans', 3 => 'Kinetoplasts', 4 => 'Other Protists'});


###############
# MANAGE output

my $log=undef;
if ($outfile) {
  open($log, '>', $outfile."_report.txt") or die "Could not open file '$outfile' $!";
}


my $log_tree=undef;
if ($outfile) {
  open(my $tree, '>', $outfile."_tree.nhx") or die "Could not open file '$outfile' $!";
  $log_tree = Bio::TreeIO->new(-fh => $tree, -format => 'nhx');
}
my $screenDisplayTree = Bio::TreeIO->new(
    -format => 'nhx',
    -fh => \*STDOUT,
    );

#############



###############
# CREATE QUERY
my $query="";
if($taxid){
	$query="txid".$taxid."[orgn]";
}
else{
	### KINGDOM LEVEL ###
	print "Please chose a Kingdom:\n";
	foreach my $key (sort{$a <=> $b} keys %KINGDOM){
		print "$key $KINGDOM{$key}\n";
	}
	print "choice:";
	my $resu_kingdom = <STDIN>;
	$resu_kingdom = key_exits($resu_kingdom, %KINGDOM);

	print "Please chose a Group within $KINGDOM{$resu_kingdom}:\n";
	foreach my $key (sort{$a <=> $b} keys %{$GROUP{$KINGDOM{$resu_kingdom} } }){
		print "$key $GROUP{$KINGDOM{$resu_kingdom}}{$key}\n";
	}
	print "choice:";
	my $resu_group=<STDIN>;
	$resu_group = key_exits($resu_group, %{$GROUP{$KINGDOM{$resu_kingdom} } });

	#case all in kingdom
	if($resu_group == "1"){
		$query=$KINGDOM{$resu_kingdom}."[Organism]";
	}
	else{ # we continue
		### GROUP LEVEL ###
		print "Please chose a Group within $GROUP{$KINGDOM{$resu_kingdom}}{$resu_group}:\n";
		foreach my $key (sort{$a <=> $b} keys %{$SUBGROUP{$GROUP{$KINGDOM{$resu_kingdom}}{$resu_group} } }){
			print "$key $SUBGROUP{$GROUP{$KINGDOM{$resu_kingdom}}{$resu_group}}{$key}\n";
		}
		print "choice:";
		my $resu_subgroup=<STDIN>;
		$resu_subgroup = key_exits($resu_subgroup, %{$SUBGROUP{$GROUP{$KINGDOM{$resu_kingdom}}{$resu_group} } });

		if($resu_subgroup == "1"){
			$query=$KINGDOM{$resu_kingdom}."[Organism] AND ".$GROUP{$KINGDOM{$resu_kingdom}}{$resu_group}."[Organism]";
		}
		else{ # we continue
			### SUBGROUP LEVEL ###
			$query=$KINGDOM{$resu_kingdom}."[Organism] AND ".$GROUP{$KINGDOM{$resu_kingdom}}{$resu_group}."[Organism] AND ".$SUBGROUP{$GROUP{$KINGDOM{$resu_kingdom}}{$resu_group}}{$resu_subgroup}."[Organism]";
		}
	}
}

############################################
# fecth the genome database using esearch
# It will give us a list of IDs
msg("### Now fetching ids into the genome database using esearch with the query: $query");
my $esearch = Bio::DB::EUtilities->new(-eutil      => 'esearch',
    		                               -email      => 'me@nbis.se',
    		                               -db         => 'genome',
                  									   -retmax 	   => [100000],
                  									   -term  => $query);
my $count = 0;
$count = $esearch->get_count;

msg("Found " . $count . " hits for " . $query . " in database '" . 'genome db' . "'\n");

if ($count == 0){ # Skip if nothing was found
	msg("Nothing found"); exit ;
}

# Go trough the XML response to extract the ids
my $xml_data;
my @id_Genomes;
$esearch->get_Response(-cb => sub { ($xml_data) = @_; } );

my $xmldoc = XML::LibXML->load_xml(string => $xml_data);

my @nodes = $xmldoc->getElementsByLocalName('Id');
foreach my $node (@nodes){
	my $avalue = $node->textContent;
	push @id_Genomes, $avalue;
}


############################################
# link genome database with taxonomy database to fetch the taxid from the id using elink
# It will give us a list of IDs
msg("### Now translating ids into taxids using elink:");
my @taxid_Genomes;
foreach my $id (@id_Genomes){
  sleep 0.5; #sleep one second. After december 2018 any site (IP address) posting more than 3 requests per second to the E-utilities without an API key will receive an error message. By including an API key, a site can post up to 10 requests per second by default... see https://www.ncbi.nlm.nih.gov/books/NBK25497/
  print "search for id = $id\n" if $verbose;
	my $elink = Bio::DB::EUtilities->new(-eutil      => 'elink',
      		                             -email      => 'me@nbis.se',
      		                             -dbfrom	   => 'genome',
      		                             -db         => 'taxonomy',
                    									 -retmax 	   => [100],
                    									 -id         => $id);

  my $xml_id = $elink->get_Response()->content();

  try{
    my $xmldoc = XML::LibXML->load_xml(string => $xml_id);
    my @nodes = $xmldoc->getElementsByLocalName('Link');
    my $taxid = undef;

    	foreach my $node (@nodes){
    		my @childnodes = $node->childNodes();

        foreach my $childnode (@childnodes){
    			my $name = $childnode->nodeName;

          if($name eq "Id"){
    				$taxid = $childnode->textContent;
    				msg("id $id has been mapped to taxid $taxid");
    				push @taxid_Genomes, $taxid;
    				last;
    			}
    		}
    		last if($taxid);
    	}
    }
    catch{
      warn "caught error: $_" if $verbose;
      msg("No taxid found for id $id");
    }
}

###############
# CREATING TREE
msg("### Now translating taxids into scientific_name using entrez:");
my $db = Bio::DB::Taxonomy->new(-source => 'entrez');
my @species_names;
my $speciesTreeReady;
foreach my $taxid (@taxid_Genomes){
    my $taxon = $db->get_taxon(-taxonid => "$taxid");
    my $spName=$taxon->scientific_name;
    msg("taxid $taxid = $spName");
    push(@species_names, $spName);
}
msg("### Now creating tree");
$speciesTreeReady = $db->get_tree(@species_names);
$speciesTreeReady->contract_linear_paths();
msg( "This is the tree of species that have whole genome sequenced:" );

# PRINT THE TREE INTO A VARIABLE
my $tree;
do {
    local *STDOUT;
    open STDOUT, ">>", \$tree;
    $screenDisplayTree->write_tree($speciesTreeReady);
};
msg($tree);
# print in a file if asked
if($log_tree){
	$log_tree->write_tree($speciesTreeReady);
}

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
  print $log $line if $log;
  print STDERR $line unless $quiet;
}

sub runcmd {
  msg("Running:", @_);
  system(@_)==0 or err("Could not run command:", @_);
}

sub key_exits{
	my ($resu_kingdom, %HASH)=@_;

	while ($resu_kingdom){
		foreach my $key (keys %HASH){
			if($resu_kingdom == $key){
				return $key;
			}
		}
		print "Wrong choice, please try again:";
		$resu_kingdom=<STDIN>;
	}


}

__END__
=head1 NAME

gaas_ncbi_get_genome_tree.pl

=head1 DESCRIPTION

The script creates a tree that covers only whole genomes from the genome NCBI database.
The result is written to the specified output file, or to STDOUT.

=head1 SYNOPSIS

    gaas_ncbi_get_genome_tree.pl [ -o outfile ]
    gaas_ncbi_get_genome_tree.pl --help

=head1 OPTIONS

=over 8

=item B<-t> or <--taxid>

To specify a specific taxid. Allow to focus on a specific part of the tree of life.

=item B<-v>

For debugging purpose.

=item B<-q>

Quiet to avoid printing the progress on STDOUT

=item B<-o>, B<--output> or B<--outfile>

The name of the output file. By default the output is the standard output

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
