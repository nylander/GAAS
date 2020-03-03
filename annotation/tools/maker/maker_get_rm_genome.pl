#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Scalar::Util qw(openhandle);
use Time::Piece;
use Time::Seconds;
use File::Basename;
use Cwd;
use GAAS::GAAS;

my $header = get_gaas_header();
my $dir = getcwd;
my $outfile = "genome.rm.fa";
my $out_fh = undef;
my $maker_dir = undef;
my $datastore = undef;
my $in = undef;
my $opt_help = 0;

# OPTION MANAGMENT
if ( !GetOptions( "help|h" => \$opt_help,
				  				"i=s" => \$in,
    			  			"outfile|o=s" => \$outfile) )
{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}
if ($opt_help) {
    pod2usage( { -verbose => 99,
                 -exitval => 0,
                 -message => "$header\n" } );
}

#######################
### MANAGE OPTIONS ####
#######################

# MANAGE IN
my @inDir;

if(! $in){
	# Find the datastore index
	my $maker_dir = undef;

	opendir(DIR, $dir) or die "couldn't open $dir: $!\n";
	my @dirList = readdir DIR;
	closedir DIR;

	my (@matchedDir) = grep $_ =~ /^.*\.maker\.output$/ , @dirList ;

	foreach my $makerDir (@matchedDir){
		push(@inDir, $makerDir);
	}
}
else{
	if (! -d "$in") {
		die "The outdirectory $in doesn't exist.\n";
	}
	else{
		push(@inDir, $in);
	}
}

# MANAGE OUT
if (-f $outfile) {
	die "The outfile $outfile already exists, exiting\n";
}
else{
	open($out_fh, '>', "$outfile") or die "Could not open file 'outfile' $!";
}

# MESSAGES
my $nbDir=$#inDir+1;
if ($nbDir == 0){die "There seems to be no maker output directory here, exiting...\n";}
print "We found $nbDir maker output directorie(s):\n";
foreach my $makerDir (@inDir){
		print "\t+$makerDir\n";
}
if ($nbDir > 1 ){print "Results will be merged together !\n";}

                #####################
                #     MAIN          #
                #####################

#############################
# Read the genome_datastore #
#############################

foreach my $makerDir (@inDir){
	my $prefix = $makerDir;
	$prefix =~ s/\.maker\.output.*//;
	my $maker_dir_path = $dir . "/" . $makerDir."/";
	my $datastore = $maker_dir_path.$prefix."_datastore" ;

	if (-d $datastore ) {
        	print "Found datastore in $makerDir, merging query.masked.fasta files now...\n";
	} else {
	        die "Could not find datastore index ($datastore), exiting...\n";
	}

	collect_recursive($datastore);

	#Close file_handler opened
	close $out_fh;
}


sub collect_recursive {
    my ($full_path) = @_;

	my ($name,$path,$suffix) = fileparse($full_path,qr/\.[^.]*/);

    if( ! -d $full_path ){

    	###################
    	# deal with fasta #
    	if($name eq "query.masked" and $suffix eq ".fasta"){

	    	#print
	    	open(my $fh, '<:encoding(UTF-8)', $full_path) or die "Could not open file '$full_path' $!";
			while (<$fh>) {
				print $out_fh $_;
			}
			close $fh;
    	}
    	return;
    }

    opendir my $dh, $full_path or die;
    while (my $sub = readdir $dh) {
        next if $sub eq '.' or $sub eq '..';

        collect_recursive("$full_path/$sub");
    }
    close $dh;
    return;
}


__END__

=head1 NAME

gaas_maker_get_rm_genome.pl

=head1 DESCRIPTION

Must be executed in the folder from which Maker was run and will find the maker output
on its own and create a concatenated masked assembly.

=head1 SYNOPSIS

    gaas_maker_get_rm_genome.pl -i maker_output_folder [-o GenomeMasked.fa]

=head1 OPTIONS

=over 8

=item B<-i>

The path to the input directory. If none given, we assume that the script is launched where Maker was run. So, in that case the script will look for the folder
*.maker.output.

=item B<--outfile>, B<-o>

The name of the masked genome file. By default, the name will genome.rm.fa

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
