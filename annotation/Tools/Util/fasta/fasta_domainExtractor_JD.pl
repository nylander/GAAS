#!/usr/bin/env perl

###################################################
# domainExtractor.pl - Jacques Dainat 01/2015     #
# Bioinformatics Infrastructure for Life Sciences #
# jacques.dainat@bils.se			  #
###################################################


use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use IO::File;

my $inputFile;
my $outputFile;
my $nameSeq;
my $start;
my $end;
my $opt_help = 0;

Getopt::Long::Configure ('bundling');
if ( !GetOptions (  'i|input_file=s' => \$inputFile,
					'n|name=s' => \$nameSeq,
				    'o|output=s' => \$outputFile,
				    's|start=i' => \$start,
					'e|end=i' => \$end,
					'h|help!'         => \$opt_help )  )
{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

if ($opt_help) {
    pod2usage( { -verbose => 2,
                 -exitval => 0 } );
}

if ((!defined($inputFile)) || (!defined($start)) || (!defined($end)) ){
	 pod2usage( { -message => 'at least 3 parameters are mandatory: -i -s and -e',
                 -verbose => 1,
                 -exitval => 1 } );
}

my $ostream     = IO::File->new();
my $ref_istream = IO::File->new();

# Manage input fasta file
$ref_istream->open( $inputFile, 'r' ) or
  croak(
     sprintf( "Can not open '%s' for reading: %s", $inputFile, $! ) );

# Manage Output
if(defined($outputFile))
{
$ostream->open( $outputFile, 'w' ) or
  croak(
     sprintf( "Can not open '%s' for reading: %s", $outputFile, $! ) );
}
else{
	$ostream->fdopen( fileno(STDOUT), 'w' ) or
      croak( sprintf( "Can not open STDOUT for writing: %s", $! ) );
}

my $nbSeq=0;
my $seq="";
my $headerCurrent="";
my $header="";
my $headerFound="no";
while (my $line = readline(*$ref_istream)) {
	$line=~ s/\s//g;
	if($line =~ m/^>/){
		$nbSeq++; $header=$headerCurrent; $headerCurrent=$line;
		if( ($nbSeq > 1) && (!defined($nameSeq)) ){
			print "The input file is an MultiFasta file.\nPlease specify the name of the sequence you are interested in (-n <name>)\n";exit;
			}
		elsif( ($nbSeq > 1) && (defined($nameSeq)) ) {
			if(($nameSeq eq $header) or (">".$nameSeq eq $header)){
				$headerFound="yes";
				last;
			}
		$seq="";
		}
	}
	else{
		if($nbSeq > 0){ #first sequence encountered. We can begin to collect sequence data
			$seq.=$line;
			}
	}
}
$header=$headerCurrent;
if($headerFound eq "no"){
	if( ($nbSeq > 1) && (!defined($nameSeq)) ){
		print "The input file is an MultiFasta file.\nPlease specify the name of the sequence you are interested in (-n <name>)\n";exit;
	}
	elsif( ($nbSeq > 1) && (defined($nameSeq)) ){
		if(($nameSeq eq $header) or (">".$nameSeq eq $header)){
			$headerFound="yes";
		}
	}
	elsif ($nbSeq == 1){
		$header=$headerCurrent;
	}
}

if( ($nbSeq > 1) && (defined($nameSeq)) && ($headerFound eq "no") ){
		print "The header you specified >$nameSeq< doesn't exist in this MultiFasta file.\nPlease check it.\n";exit;
} 

print "Name studied sequence: $header\n";
print "sequence: $seq\n";
if($start<0 || $end <0){print "Start and End cannot be a negative value!\n"; exit;}
if(length($seq) < $start){print "Start position for extraction is over the sequence size !\n"; exit;}
if(length($seq) < $end){print "End position for extraction is over the sequence size !\n"; exit;}
#end is 1-based coordinate system and 0-based coordinate system
#start is 1-based coordinate system
# The extraction compute in 0-based coordinate system
# Lets change the 1-based coordinate system in 0-based coordinate system for the start
$start=$start-1;
my $lengtExtraction=$end-$start; #Length in 0-based coordinate (in 1-based coordinate we must add +1)
print "Length sequence extracted: $lengtExtraction\n";
my $extractedPart=substr($seq, $start, $lengtExtraction);
print "Sequence extracted: $extractedPart\n";


__END__

=head1 NAME

domainExtractor_JD.pl -
The script allows to extract a part of a AA or nt sequence.
The script takes as input a (multi)fasta file and coordinates of part that you want extract.
If the Input file is a MultiFastaFile you have to specify to the script the header of the sequence you want to extract.
NOTE: The script expect the use of 1-based coordinate system. So, -s 1 -e 1 extract the first AA/nt
/!\ Some file formats are 1-based (GFF, SAM, VCF) and others are 0-based (BED, BAM)
/!\ Ensembl uses 1-based coordinate system when UCSC uses 0-based coordinate system
/!\ Be aware of what kind of coordinate you are using as input.
Rule of coordinate system
	1-based coordinate system = Numbers nucleotides directly
    0-based coordinate system = Numbers between nucleotides
    
=head1 SYNOPSIS

    domainExtractor_JD.pl -i <input file> -s <start_coordinate> -e <end_coordinate> [-o <output file> -n <sequence name>]
    domainExtractor_JD.pl --help

=head1 OPTIONS

=over 8

=item B<-i>, B<--file> or B<-ref>

Input fasta file that will be read.

=item B<-s> or B<--start> 

Start coordinate of the region that will be extract

=item B<-e> or B<--end> 

End coordinate of the region that will be extract

=item B<-n> or B<--name> 

In Multifasta file case, the name allows to specify which sequence you are interested in.

=item B<-o> or B<--output> 

Output file.  If no output file is specified, the output will be
written to STDOUT.



=back

=cut
