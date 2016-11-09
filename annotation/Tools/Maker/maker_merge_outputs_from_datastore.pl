#!/usr/bin/env perl

use strict;
use File::Copy;
use warnings;
use Scalar::Util qw(openhandle);
use Time::Piece;
use Time::Seconds;
use Cwd;
use Pod::Usage;
use Getopt::Long qw(:config no_ignore_case bundling);
use BILS::Handler::GXFhandler qw(:Ok);
use BILS::Handler::GFF3handler qw(:Ok);
use Bio::Tools::GFF;
use IO::File;
use File::Basename;

my $header = qq{
########################################################
# BILS 2015 - Sweden                                   #  
# jacques.dainat\@bils.se                               #
# Please cite BILS (www.bils.se) when using this tool. #
########################################################
};

my $out = undef;
my $in = undef;
my $help= 0;

if ( !GetOptions(
    "help|h" => \$help,
    "i=s" => \$in,
    "output|outfile|out|o=s" => \$out))

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
if(! $out){
	$out="annotations";
}

# STANDARD
my $protein_file = "annotations.proteins.fa";
my $annotations_file = "annotations.gff";

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
if (-d "$out") {
	print "The output directory <$out> already exists, we skip the merge step.\n";
} 
else{
	mkdir $out;
	open(my $gff_out, '>', $out."/".$annotations_file) or die "Could not open file '$out/$annotations_file' $!";
	open(my $protein_out, '>', $out."/".$protein_file) or die "Could not open file '$out/$protein_file' $!";

	foreach my $makerDir (@inDir){
		my $prefix = $makerDir;
		$prefix =~ s/\.maker\.output//;
		my $maker_dir_path = $dir . "/" . $makerDir."/";
		my $datastore = $maker_dir_path.$prefix."_datastore" ;

		if (-d $datastore ) {
	        	print "Found datastore in $makerDir, merging annotations now...\n";
		} else {
		        die "Could not find datastore index ($datastore), exiting...\n";
		}

		# This is one way to open a file...
		#open (my $IN, '<', $datastore) or die "FATAL: Can't open file: $datastore for reading.\n$!\n";

		#Read list of dir
		opendir(DIR, $datastore) or die "cannot open dir $datastore: $!";
	  	my @listDirL1= readdir DIR;
	  	closedir DIR;

	  	foreach my $dirL1 (@listDirL1){
	  		if($dirL1 =~ /^\./){next;}

	  		#Read list of dir
			opendir(DIR, $datastore."/".$dirL1) or die "cannot open dir $datastore/$dirL1: $!";
	  		my @listDirL2= readdir DIR;
	  		closedir DIR;

	  		foreach my $dirL2 (@listDirL2){
	  			if($dirL2 =~ /^\./){next;}

		  		#Read list of dir
				opendir(DIR, $datastore."/".$dirL1."/".$dirL2) or die "cannot open dir $datastore/$dirL1/$dirL2: $!";
		  		my @listDirL3= readdir DIR;
		  		closedir DIR;

		  		foreach my $dirL3 (@listDirL3){
		  			if($dirL3 =~ /^\./){next;}

		  			opendir(DIR, $datastore."/".$dirL1."/".$dirL2."/".$dirL3) or die "cannot open dir $datastore/$dirL1/$dirL2/$dirL3: $!";
		  			my @listFile= readdir DIR;
		  			closedir DIR;

		  			my $gff_file = undef;
		  			my $aa_file = undef;
		  			if(grep $_ =~ /\.gff$/ , @listFile ){
		  				my @matchFile = grep $_ =~  /\.gff$/, @listFile ;
						$gff_file = $matchFile[0] ;
		  			}

		  			if(grep $_ =~ /\.proteins\.fasta$/ , @listFile ){
		  				my @matchFile = grep $_ =~  /\.proteins\.fasta$/, @listFile ;
						$aa_file = $matchFile[0] ;
		  			}

					if ($gff_file) {
						my $filename="$datastore/$dirL1/$dirL2/$dirL3/$gff_file";
						open(my $fh, '<:encoding(UTF-8)', $filename) or die "Could not open file '$filename' $!";
						while (<$fh>) {
							print $gff_out $_;
						}
						close $fh;
					}
				
					if ($aa_file ) {
						my $filename="$datastore/$dirL1/$dirL2/$dirL3/$aa_file";
						open(my $fh, '<:encoding(UTF-8)', $filename) or die "Could not open file '$filename' $!";
						while (<$fh>) {
							print $protein_out $_;
						}
						close $fh;
					}
				}
			}
		}
	}

	close $gff_out;
	close $protein_out;
}

print "Now save a copy of the Maker option files ...\n";
if (-f "$out/maker_opts.ctl") {
	print "A copy of the Maker files already exists in $out/maker_opts.ctl.  We skip it.\n";
}
else{
	if(! $in){
		copy("maker_opts.ctl","$out/maker_opts.ctl") or print "Copy failed: $! \n";
		copy("maker_exe.ctl","$out/maker_exe.ctl") or print "Copy failed: $! \n";
		copy("maker_evm.ctl","$out/maker_evm.ctl") or print "Copy failed: $! \n";
		copy("maker_bopts.ctl","$out/maker_bopts.ctl") or print "Copy failed: $! \n";
	}
	else{
		my ($name,$path,$suffix) = fileparse($in);
		copy("$path/maker_opts.ctl","$out/maker_opts.ctl") or print  "Copy failed: $!";
		copy("$path/maker_exe.ctl","$out/maker_exe.ctl") or print "Copy failed: $!";
		copy("$path/maker_evm.ctl","$out/maker_evm.ctl") or print "Copy failed: $!";
		copy("$path/maker_bopts.ctl","$out/maker_bopts.ctl") or print "Copy failed: $!";
	}
}


############################################
# Now manage to split file by kind of data # Split is done on the fly (no data saved in memory)
############################################
my $splitedData_dir= "$out/annotationByType";

if (-f $splitedData_dir) {
	print "Output $splitedData_dir of the <split by type> step already exists. We skip it."
}
else{
	print "Now split file by data type...\n";
	mkdir $splitedData_dir;
	# split data by column 2 with awk
	exec "awk '{if(\$2 ~ /[a-zA-Z]+/) print \$0 > \"$splitedData_dir/\"\$2\".gff\"}' $out\"/\"$annotations_file";
}

print "All done!\n";

# --------------



=head1 NAME

maker_merge_outputs.pl -
Usage:
	Must be executed in the folder from which Maker was run and will find the maker output
	on its own and create a concatenated annotation file. 

=head1 SYNOPSIS

    ./maker_merge_outputs.pl 
    ./maker_merge_outputs.pl --help

=head1 OPTIONS

=over 8

=item B<-i>

The path to the input directory. If none given, we assume that the script is launched where Maker was run. So, in that case the script will look for the folder 
*.maker.output.

=item B<-o> or B<--output>

The name of the output directory. By default the name is annotations

=item B<-h> or B<--help>

Display this helpful text.

=back

=cut
