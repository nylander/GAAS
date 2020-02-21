#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Pod::Usage;
use Scalar::Util qw(openhandle);
use Time::Piece;
use Time::Seconds;
use File::Basename;
use Cwd;

my $header = qq{
########################################################
# NBIS 2018 - Sweden                                   #  
# jacques.dainat\@nbis.se                               #
# Please cite NBIS (www.nbis.se) when using this tool. #
########################################################
};

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
    pod2usage( { -verbose => 2,
                 -exitval => 2,
                 -message => "$header\n" } );
}


#######################
### MANAGE OPTIONS ####
#######################


# MANAGE IN
my @inDir;
my $dir = getcwd;

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


# Streaming the file, line by line
while (<$IN>) {
	chomp; # Trims the line, removes the line breaks
	my $line = $_; # store the line in a variable
	next unless ($line =~ /^.*FINISHED$/) ; # We only want finished contig annotations...
	
	my ($contig,$location,$status) = split("\t",$line);
	
	# If the contig includes a dot character, the output file will include a percent character...
	if ($contig =~ /^.*\..*$/) {
		$contig =~ s/\./\%2E/g ;
	}

	my $void = $location . "theVoid." . $contig ;
	
	my $repeat_contig = $maker_dir . "/" . $void . "/" . "query.masked.fasta" ;
	
	if (-f $repeat_contig ) {
		system("cat $repeat_contig >> $outfile");
	}
	
}
# We should close the file to make sure that the transaction finishes cleanly.
close ($IN);


=head1 NAME

maker_get_rm_genome.pl

Must be executed in the folder from which Maker was run and will find the maker output
on its own and create a concatenated masked assembly.

=head1 SYNOPSIS

    ./maker_get_rm_genome.pl

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

=cut
