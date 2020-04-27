#!/usr/bin/env perl


#
# NCBI recommends that users post no more than three URL requests per second
#

use strict;
use Try::Tiny;
use Getopt::Long;
use Bio::DB::EUtilities;
use Bio::DB::Fasta;
use Bio::SeqIO;
use Scalar::Util qw(openhandle);
use Time::Piece;
use Time::Seconds;
use Pod::Usage;
use XML::LibXML;
use GAAS::GAAS;

my $header = get_gaas_header();
my $opt_output = undef;
my $col = undef;
my $message="";
my $list=undef;
my $separator=undef;
my $lineToAvoid=undef;
my $help;
my $quiet = undef;

if ( !GetOptions(
    "help|h" 							=> \$help,
		"list|l=s" 						=> \$list,
		"line=i" 							=> \$lineToAvoid,
		"col=i" 							=> \$col,
		"s=s"									=>\$separator,
		"q" 									=> \$quiet,
    "o|output|outfile=s" => \$opt_output))
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

if ( ! (defined($list)) ){
    pod2usage( {
           -message => "$header\nAt least 1 parameters is mandatory (--list).\n",
           -verbose => 0,
           -exitval => 1 } );
}
###############
# MANAGE output

my $log=undef;
if ($opt_output) {
  open($log, '>', $opt_output."_log.txt") or die "Could not open file '$opt_output'_log.txt $!";
}
my $error=undef;
if ($opt_output) {
  open($error, '>', $opt_output."_error.txt") or die "Could not open file '$opt_output'_error.txt $!";
}

my $outstream;
if ($opt_output) {
  $opt_output=~ s/.fasta//g;
  $opt_output=~ s/.fa//g;
  open($outstream, '>', $opt_output.".fa") or die "Could not open file '$opt_output' $!";
 # $ostream= Bio::SeqIO->new(-fh => $fh, -format => 'Fasta' );
}
else{
  $outstream = \*STDOUT;
}

#Manage column with the ID
if (! defined $col){
	$col=0;
}
else{$col=$col -1 ;}

#Manage line to avoid
if (! defined $lineToAvoid){
	$lineToAvoid=0;
}

print "The ID will be retrieved into the column ".($col+1)." of the file $list\n";
print "The first $lineToAvoid lines will be ignored.\n";

########################
### Manage INPUT FILES #
my %list_of_ID;
my $nbID=0;
if (-f $list){
	my $ID_list  = IO::File->new("<".$list);

	##############
	#### MAIN ####

	# create hash of ID to remove
	my $cpt_line=0;
	while ( <$ID_list> ) {
	  $cpt_line++;

	  if($cpt_line > $lineToAvoid){

		  chomp;
		  if(! $_ =~ /^\s*$/){

			my @cols;
		  	if (! $separator){
		  		@cols = split /\s/, $_;
		 	}
		 	else{
		 		@cols = split /$separator/, $_;
		 	}

	        my $id = $cols[$col];
	        $id =~ s/[^[:print:]]+//g;
	        print $id."\n";
	        $list_of_ID{$id}++;
	        $nbID++;
		  }
		}
	}
}

msg( "There is $nbID identifiers, lets try to extract the corresponding sequences\n");

##Fetch correct ID
foreach my $ID (keys %list_of_ID){

	############################################
	# fecth the genome database using esearch
	# It will give us a list of IDs
	msg("### Now fetching the correct identifier into the protein database using esearch with the query: $ID\n");
	my $idcorrect=undef;
	try{
		my $factory = Bio::DB::EUtilities->new(-eutil      => 'esearch',
				                               -email      => 'me@foo.com',
				                               -db         => 'protein',
											   -retmax 	   => [100000],
											   -term  => $ID);
		my $count = 0;
		$count = $factory->get_count;

		if ($count == 0){ # Skip if nothing was found
			if ($opt_output) {
				print $error "a - No identifier found for $ID\n"; ## => print to log error
			}
			else {
				msg("a - No identifier found for $ID\n");
			}
			next;
		}
		else{
			#msg("We found $count ID for $ID  in database 'protein db' \n");
		}

		# Go trough the XML response to extract the ids
		my $xml_data;
		$factory->get_Response(-cb => sub { ($xml_data) = @_; } );

		my $xmldoc = XML::LibXML->load_xml(string => $xml_data);

		my @nodes = $xmldoc->getElementsByLocalName('Id');


		foreach my $node (@nodes){
			$idcorrect = $node->textContent;
			last;
		}
		sleep(1)
	}
	catch {
		warn "caught error: $_"; # not $@
	};

	if(! $idcorrect){
		msg("b - No identifier found for $ID\n");
		print $error "b - No identifier found for $ID\n"; ## => print to log error
	}
	else{
	##Fetch sequence from correct ID

		############################################
		# fecth the genome database using esearch
		# It will give us a list of IDs
		msg("### Now fetching the sequence into the protein database using efetch with the query: $idcorrect\n");
		my $factory = Bio::DB::EUtilities->new(-eutil      => 'efetch',
				                               -email      => 'me@foo.com',
				                               -db         => 'protein',
											   -retmax 	   => [100000],
											   -rettype => 'fasta',
											   -id  => $idcorrect);

		my $fasta =	$factory->get_Response->content;
		print $outstream $fasta;
	}
}

if ($opt_output) {
	close $log;
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


=head1 NAME

gaas_ncbi_get_sequence_from_list.pl

=head1 DESCRIPTION


The script allow to retrieve the sequences from the NCBI ID list.
The list should be a column in a file containing one ID per line.
The result is written to the specified output file, or to STDOUT in fasta format.

=head1 SYNOPSIS

    gaas_ncbi_get_sequence_from_list.pl --list file.txt [ -o outfile ]
    gaas_ncbi_get_sequence_from_list.pl --help

=head1 OPTIONS

=over 8

=item B<--list> or or B<-l>

File containing ID by colomn

=item B<--line>

Integer, number of line to avoid. Allow to avoid headers.

=item B<--col>

column containing the ID. By default the first column is considered.

=item B<-q>

Field separator, by default un-printable character are use as separator (\s). You can define the one you wnat with this option.

=item B<-q>

Quiet to avoid any print on STDOUT

=item B<-o>, B<--output> or B<--outfile>

The name of the output file. By default the output is the standard output

=item B<-h> or B<--help>

Display this helpful text.

=back

=cut
