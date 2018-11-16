#!/usr/bin/env perl

use strict;
use File::Copy;
use warnings;
use Scalar::Util qw(openhandle);
use Time::Piece;
use Time::Seconds;
use Cwd;
use Pod::Usage;
use URI::Escape;
use Getopt::Long qw(:config no_ignore_case bundling);
use BILS::Handler::GXFhandler qw(:Ok);
use BILS::Handler::GFF3handler qw(:Ok);
use Bio::Tools::GFF;
use IO::File;
use File::Basename;
use IPC::Cmd qw[can_run run];

my $header = qq{
########################################################
# NBIS 2015 - Sweden                                   #  
# jacques.dainat\@nbis.se                               #
# Please cite NBIS (www.nbis.se) when using this tool. #
########################################################
};

my %file_hds;
my $out = undef;
my $in = undef;
my $help= 0;

if ( !GetOptions(
    "help|h" => \$help,
    "i=s" => \$in,
    "output|out|o=s" => \$out))

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

# STANDARD
my $default_out_dir_name = "maker_output_processed";
my $outfile = "genome";

# MANAGE OUT
if(! $out){
	$out = $default_out_dir_name;
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
if (-d "$out") {
	print "The output directory <$out> already exists, we skip the merge step.\n";
} 
else{
	print "Creating the $out folder\n";
	mkdir $out;

	foreach my $makerDir (@inDir){
		my $prefix = $makerDir;
		$prefix =~ s/\.maker\.output.*//;
		my $maker_dir_path = $dir . "/" . $makerDir."/";
		my $datastore = $maker_dir_path.$prefix."_datastore" ;

		if (-d $datastore ) {
	        	print "Found datastore in $makerDir, merging annotations now...\n";
		} else {
		        die "Could not find datastore index ($datastore), exiting...\n";
		}

		if ( grep -f, glob '$out/*.fasta' ) {
			print "Output fasta file already exists. We skip the fasta gathering step.\n";
		}
		elsif ( grep -f, glob '$out/*.gff' ) {
			print "Output gff file already exists. We skip the gff gathering step.\n";
		}
		else{
			collect_recursive($datastore);

			#Close all file_handler opened
			foreach my $key (keys %file_hds){
				close $file_hds{$key};
			}
		}
	}
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

#make the annotation safe
my $annotation="$out/maker.gff";
if (-f $annotation) {
	print "Protecting the maker.gff annotation by making it readable only\n";
	system "chmod 444 $annotation";
}
else{
	print "ERROR: Do not find the $annotation file !\n";
}


#do statistics
my $full_path = can_run('gff3_sp_statistics.pl') or print "Cannot launch statistics. gff3_sp_statistics.pl script not available\n";
if ($full_path) {
        print "Performing the statistics of the maker.gff annotation file\n";
	my $annotation_stat="$out/maker_stat.txt";
        system "gff3_sp_statistics.pl --gff $annotation -o $annotation_stat";
}

print "All done!\n";

sub collect_recursive {
    my ($full_path) = @_;
	
	my ($name,$path,$suffix) = fileparse($full_path,qr/\.[^.]*/);

    if( ! -d $full_path ){
    	
    	###################
    	# deal with fasta #
    	if($suffix eq ".fasta"){
    		my $key = undef;
    		my $type = undef;
    		if($name =~ /([^\.]+)\.transcripts/){
    			$key = $1;
    			$type = "transcripts";
    		}
    		if($name =~ /([^\.]+)\.proteins/){
    			$key = $1;
    			$type = "proteins";
    		}
    		if($name =~ /([^\.]+)\.noncoding/){
    			$key = $1;
    			$type = "noncoding";
    		}
    		if($key){
    			my $source = 'maker';
	    		$source .= ".$key" if($key ne 'maker');
	    		
				my $prot_out_file_name = "$outfile.all.$source.$type.fasta";
				my $protein_out_fh=undef;
				if( _exists_keys (\%file_hds,($prot_out_file_name)) ){
					$protein_out_fh = $file_hds{$prot_out_file_name};
				}
				else{
					open($protein_out_fh, '>', "$out/$prot_out_file_name") or die "Could not open file '$out/$prot_out_file_name' $!";
					$file_hds{$prot_out_file_name}=$protein_out_fh;
				}
				 	
	    		#print
	    		open(my $fh, '<:encoding(UTF-8)', $full_path) or die "Could not open file '$full_path' $!";
					while (<$fh>) {
						print $protein_out_fh $_;
					}
				close $fh;
    		}
    	}

    	################
 		#deal with gff #
    	if($suffix eq ".gff"){
			system "awk '{if(\$2 ~ /[a-zA-Z]+/) { gsub(/:/, \"_\" ,\$2); print \$0 >> \"$out/\"\$2\".gff\"}}' $full_path";
    	}
    	
    	return;
    }
    if($name =~ /^theVoid/){ # In the void there is sub results already stored in the up folder. No need to go such deep otherwise we will have duplicates.
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

sub _exists_keys {
    my ($hash, $key, @keys) = @_;

    if (ref $hash eq 'HASH' && exists $hash->{$key}) {
        if (@keys) {
            return exists_keys($hash->{$key}, @keys);
        }
        return 1;
    }
    return '';
}

__END__
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
