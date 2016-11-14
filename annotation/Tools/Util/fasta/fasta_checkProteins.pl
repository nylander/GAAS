#!/usr/bin/env perl

#############################################
# checkProtein.pl - Jacques Dainat 12/2014 # 
#############################################

use strict;
use warnings;
use Getopt::Long;
use IO::File;
use Pod::Usage;

#VERIABLE DECLARATION
my $opt_reffile;
my $opt_help;
my $opt_output;
my $nbProt=0;

# OPTION MANAGMENT
if ( !GetOptions( 'f|ref|reffile=s' => \$opt_reffile,
				  'o|output=s'      => \$opt_output,
				  'h|help!'         => \$opt_help ) )

{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

if ($opt_help) {
    pod2usage( { -verbose => 2,
                 -exitval => 0 } );
}
 
if ( ! (defined($opt_reffile)) ){
    pod2usage( {
           -message => "\nAt least 1 parameter is mandatory:\nInput fasta file (--f)\n\n",
           -verbose => 0,
           -exitval => 2 } );
}

### MANAGE OUTPUT
my $ostream     = IO::File->new();
if ($opt_output){
        $ostream->open( $opt_output, 'w' ) or
        croak(
            sprintf( "Can not open '%s' for writing %s", $opt_output."/GOFeatures.gff", $! )
        );
}
else{
	
	   $ostream->fdopen( fileno(STDOUT), 'w' ) or
	        croak( sprintf( "Can not open STDOUT for writing: %s", $! ) );
	 
}
my $output= $ostream ;

### READ FILE
if (! (-f $opt_reffile)){
	print "File doesnot exist\n"; exit;
}

open(FIC,$opt_reffile) or die "Couldn't open the file $opt_reffile\n";

my $header; my $sequence; my $firstHeaderRead="no";
my $nbProtWithStop=0; my $nbProtWithoutStop=0; my $nbProtWithStart=0; my $nbProtWithoutStart=0; my $nbProtWithStopStart=0; my $nbProtWithoutStartStop=0; my $specialStart=0;
my $nbProtWithStartWithoutStop=0;my $nbProtWithoutStartWithStop=0;
while( my $line = <FIC> ) {
        chomp($line) ;
        if($line =~ m/^>/){
            $nbProt++;
        	if ($firstHeaderRead eq "yes"){ # Allow to avoid if comment line at the beginning of the file
        		if ($sequence =~ m/^L/){
                    $specialStart++;
                }
                
        		if ($sequence =~ m/^M/){
					$nbProtWithStart++;
					if ( $sequence =~ m/[\.X\*]$/ ){
						$nbProtWithStop++;
						$nbProtWithStopStart++;
        			}
        			else{$nbProtWithStartWithoutStop++;$nbProtWithoutStop++;}

        		}
        		else{$nbProtWithoutStart++;
        			if ( $sequence =~ m/[\.X\*]$/) {
        				$nbProtWithStop++;
						$nbProtWithoutStartWithStop++;
					}
					else{$nbProtWithoutStartStop++;$nbProtWithoutStop++;}
        		}

        		$sequence="";
        	}
        	$header=$line;
        	$firstHeaderRead="yes";
        }
        elsif($firstHeaderRead eq "yes"){ 
        	$sequence.=$line;

        }
}
# Check last protein read
if ($sequence =~ m/^L/){
    $specialStart++;
}
if ($sequence =~ m/^M/){
    $nbProtWithStart++;
    if ( $sequence =~ m/[\.X\*]$/ ) {
        $nbProtWithStop++;
        $nbProtWithStopStart++;
    }
    else{$nbProtWithStartWithoutStop++;$nbProtWithoutStop++;}
}
else{$nbProtWithoutStart++;
    if ( $sequence =~ m/[\.X\*]$/ ) {
        $nbProtWithStop++;
        $nbProtWithoutStartWithStop++;
    }
    else{$nbProtWithoutStartStop++;$nbProtWithoutStop++;}
}

my $Result;
$Result = "\nWe checked $nbProt Proteins:\n";
$Result .= "M....?  We have $nbProtWithStart proteins with a start at the first position \n";
$Result .= "?....X  We have $nbProtWithStop proteins with a stop at the last position \n";
$Result .= "M....X  We have $nbProtWithStopStart proteins with a start at the first position and a stop at the last position (Correct Proteins !!)\n\n";
$Result .= ".....?  We have $nbProtWithoutStart proteins without a start at the first position \n";
$Result .= "?.....  We have $nbProtWithoutStop proteins without a stop at the last position \n";
$Result .= "......  We have $nbProtWithoutStartStop proteins without a start at the first position and without stop at the last position\n\n";
$Result .= "M.....  start wihtout stop= $nbProtWithStartWithoutStop\n";
$Result .= ".....X  stop wihtout start= $nbProtWithoutStartWithStop\n";
print "$Result" if($output);
print $output $Result; 

my $prop = ($specialStart*100)/$nbProtWithoutStart;
print "special start (L) $specialStart corresponding to $prop% of sequence without start. If this value is close to 10% it should correspond to the number we can find randomly (Mean that prediction didn't take this potential start codon in account.)\n";
__END__

=head1 NAME

checkProteins.pl -
The script take a fasta file as input. -
It will check the presence of Start (M in first position) and Stop (. or X or * at last position) of each sequence.

=head1 SYNOPSIS

    ./checkProteins.pl -f=infile.fa [ -o outfile ]
    ./checkProteins.pl --help

=head1 OPTIONS

=over 8

=item B<-f>, B<--reffile> or B<-ref>

Input fasta file that will be read. In general come from gffread output.

=item B<-o> or B<--output> 

By default the result is written on screen at te fly. If you give an output it will writte the report in this file.

=back

=cut