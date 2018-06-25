#!/usr/bin/env perl

use strict;
use warnings;
use IO::File;
use File::Basename;
use Getopt::Long;
use Pod::Usage;

my @copyARGV=@ARGV;

my $opt_output = undef;
my $datfile = undef;
my $emblfile = undef;
my $help= undef;

my $header = qq{
########################################################
# BILS 2016 - Sweden                                   #  
# jacques.dainat\@bils.se                               #
# Please cite BILS (www.bils.se) when using this tool. #
########################################################
};

if ( !GetOptions("dat=s"		  => \$datfile,
				"embl=s"		  => \$emblfile,
		    	"o|out=s" 	  => \$opt_output,
		    	"h|help"	  => \$help) )
{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

if ($help) {
    pod2usage( { -verbose => 2,
                 -exitval => 2,
                 -message => "$header\n" } );
}

if ( ! defined( $datfile) or ! defined( $emblfile)) {
    pod2usage( {
           -message => "$header\nMust specify at least 2 files. dat and embl\n",
           -verbose => 0,
           -exitval => 1 } );
}

##### Stream in 1 
my $fh1;
if ($datfile) {
  open($fh1, '<', $datfile) or die "Could not open file '$datfile' $!";
}

##### Stream in 2
my $fh2;
if ($emblfile) {
  open($fh2, '<', $emblfile) or die "Could not open file '$emblfile' $!";
}


##### Stream out
my $emblout;
if ($opt_output) {
  open($emblout, '>', $opt_output) or die "Could not open file '$opt_output' $!";
}

my %sizes;
my %headers;
my $ID = undef;
my $sourceSeen=undef;
while( my $line = <$fh1>)  {   
   
    if( $line =~ m/^ID/){
    	#
    	my @list = split(/\s/,$line);
    	$ID = $list[11];
		$sizes{$ID}++;
		$sourceSeen=undef;
    }
    
    ################## Look for signal to stop save the information #############################
    # #we keep all until source and the few next lines related to source. (Stop at another FT than source or XX if no other source available.)
    #
    if( $line =~ m/^FT   source/){
    	$sourceSeen=1;
    }

    if( $line =~ m/^FT   [^source|^\s]/){ 
    	$ID=undef;
    }
 	
 	if( $sourceSeen){
 		 if( $line =~ m/^XX/){
 		 	$ID=undef;
 		 }
 	}
	###############################################


    if($ID){
    	$headers{$ID}.=$line;
    }
}

my $uniq= 1;
foreach my $key (keys %sizes){
	my $nb = $sizes{$key};
	if ($nb != 1){
		print $key." ".$nb."\n";
		$uniq = undef;
	}
}
if ($uniq){
 print "Fine, all contig size are uniq we can use them to map the two files information.\n";
}

#everything needed from dat file is saved in headers now
#foreach my $key (keys %headers){
#		print $headers{$key};
#}

# print $fh2 but part (ID line until source and few next lines) are replaced by those saved from the dat file. (We do that if an comaon identifier is found: here the size in bp fron the ID line is the ioentifier) 
my $printNext=1;
while( my $line = <$fh2>)  {   

	if( $line =~ m/^ID/){
    	#
    	my @list = split(/\s/,$line);
    	$ID = $list[11];
		$sourceSeen=undef;
    

    	if ( exists($headers{$ID}) ){
    		print $emblout $headers{$ID};
    		$printNext=undef;
    	}
    }

    if($printNext){
    	print $emblout $line;
    }
    else{
    	if( $line =~ m/^FT   source/){
    		$sourceSeen=1;
    	}

    	if( $line =~ m/^FT   [^source|^\s]/){ #we keep all wat is related tosource as well and then we stop.
    		$printNext=1;
    	}
 	
 		if( $sourceSeen){
 			if( $line =~ m/^XX/){
 		 		$printNext=1;
 		 	}
 		}

 		if($printNext){
    		print $emblout $line;
    	}
    }

}

__END__


=head1 NAME
 
sync_dat_and_embl.pl - This script allow to update the record "headers" of an EMBL file by those from a dat file provided by ENA.
It is useful when an assembly/annotation has been submitted using AGP file while the annotation has been done directly on the chromosomes. 
(Passing by an AGP file happens only if the chromosome and unplaced (not related at all to any chromosome) contigs are part of the same assembly).


=head1 SYNOPSIS

    ./sync_dat_and_embl.pl --dat=infile --embl=infile2 -o=outFile 
    ./sync_dat_and_embl.pl --help

=head1 OPTIONS

=over 8

=item B<--dat>

Input dat file provided by ENA

=item B<--embl>

Input embl file

=item  B<--out>, B<--output> or B<-o>

The output will be the EMBL file with the record "headers" modified

=item B<--help> or B<-h>

Display this helpful text.

=back

=cut
