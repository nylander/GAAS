#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Bio::DB::EUtilities;
use Bio::SeqIO;
use Scalar::Util qw(openhandle);
use Time::Piece;
use Time::Seconds;
use Pod::Usage;
use XML::LibXML;
use Data::Dumper;

use Bio::DB::Taxonomy;
use Bio::TreeIO;
use Bio::Tree::Tree;
use Bio::Tree::TreeFunctionsI;

my $header = qq{
########################################################
# BILS 2015 - Sweden                                   #  
# Marc Hoeppner / Jacques Dainat                        #
# Please cite BILS (www.bils.se) when using this tool. #
########################################################
};

my $outfile = undef;
my $format = "fasta";
my $quiet;
my $organisms = undef;
my $dbs = undef;
my $outdir = "tmp";
my $taxid=undef;
my $list;
my $help;

if ( !GetOptions(
    "help|h" => \$help,
	"t|taxid=i" => \$taxid,
	"outdir=s" => \$outdir,
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


# .. Create output directory

runcmd("mkdir -p $outdir");

# .. set up log file

my $logfile = "$outdir/reference_sequences.log";
msg("Writing log to: $logfile");
open LOG, '>', $logfile or err("Can't open logfile");


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
	
my @taxid_Genomes;
## fecth database
my $factory = Bio::DB::EUtilities->new(-eutil      => 'esearch',
		                               -email      => 'me@foo.com',
		                               -db         => 'genome',
									   -retmax 	   => [100000],
									   -term  => $query);
my $count = 0;
$count = $factory->get_count;

msg("Found " . $count . " hits for " . $query . " in database '" . 'genome db' . "'\n"); 
		
if ($count == 0){ # Skip if nothing was found
	print "Nothing found"; exit ;
}

my $xml_data;
$factory->get_Response(-cb => sub { ($xml_data) = @_; } );

my $xmldoc = XML::LibXML->load_xml(string => $xml_data);

my @nodes = $xmldoc->getElementsByLocalName('Id');
foreach my $node (@nodes){
	my $avalue = $node->textContent;
	push @taxid_Genomes, $avalue;
}

###############
# CREATING TREE
my $db = Bio::DB::Taxonomy->new(-source => 'entrez');
my @species_names;
my $speciesTreeReady;
foreach my $taxid (@taxid_Genomes){
	print "ll $taxid\n";
        my $taxon = $db->get_taxon(-taxonid => "$taxid");
        print "$taxon\n";
        my $spName=$taxon->scientific_name;
#       print "$taxid = $spName\n";
        push(@species_names, $spName);
    }
    $speciesTreeReady = $db->get_tree(@species_names);
#    print $speciesTreeReady, "\n"; 
    ## Clean Tree
    $speciesTreeReady->contract_linear_paths();


exit;





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
