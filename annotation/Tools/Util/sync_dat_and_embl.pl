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
my $bac = undef; #believe in AC
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
                "bac!"     => \$bac,
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
my $ACline=0;
while( my $line = <$fh1>)  {   
   
    if( $line =~ m/^ID/){
    	#
    	my @list = split(/\s/,$line);
    	$ID = $list[11];
        #$line=$list[0]." ".$list[1]." ".$list[2]." ".$list[3]." ".$list[4]." ".$list[5]." ".$list[6]." ".$list[7]." ".$list[8]." STD; ".$list[10]." ".$list[11]." ".$list[12]."\n";
        #print $line."\n";
		$sizes{$ID}++;
        #print $ID."\n";
		$sourceSeen=undef;
    }
 
    ################## Look for signal to stop save the information #############################
    # #we keep all until source and the few next lines related to source. (Stop at another FT than source or XX if no other source available.)
    #
    if($bac){
        if( $line =~ m/^AC/){
            $ACline++;
            if($ACline % 2 == 0){
                my @list = split(/\s/,$line);
                my $newID= $list[2];
                chomp($newID);
                $sizes{$newID} = delete $sizes{$ID};
                $headers{$newID} = delete $headers{$ID};
                $ID=$newID;
                #print "$ID\n";
            }
        }
    }
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
        print "/!\ This size is not uniq:\n";
		print $key." ".$nb."\n";
		$uniq = undef;
	}
}
if ($uniq){
 my $nbID = keys %sizes;
 print "Fine, all contig size are uniq we can use them to map the two files information.\nThere is $nbID keys in the dat file.\n";
}

#everything needed from dat file is saved in headers now
#foreach my $key (keys %headers){
#		print $headers{$key};
#}

# print $fh2 but part (ID line until source and few next lines) are replaced by those saved from the dat file. (We do that if an comon identifier is found: here the size in bp fron the ID line is the ioentifier) 
my $printNext=1;
my $nbIDfound=0;
my $nbIDTotal=0;
$ACline=0;
my $saved_line=undef;
while( my $line = <$fh2>)  {   
    my $skip_line=undef;

	if($saved_line){
        $saved_line.=$line;
    }

    if( $line =~ m/^ID/){
    	if($bac){$printNext=undef;};
        $saved_line.=$line;

        #
    	my @list = split(/\s/,$line);
    	$ID = $list[10];
        #print "ID= $ID\n";
		$sourceSeen=undef;
    
    	if ( exists($headers{$ID}) and ! $bac){
            $nbIDfound++;
    		print $emblout $headers{$ID};
    		$printNext=undef;
    	}
    }

    if($bac){
       if( $line =~ m/^AC/){
            $ACline++;
            if($ACline % 2 == 0){
                $nbIDTotal++;
                my @list = split(/\s/,$line);
                $ID=$list[2];
                chomp($ID);
                if ( exists($headers{$ID}) ){
                    $nbIDfound++;
                    print "$ID\n";
                    print $emblout $headers{$ID};
                    $saved_line="";
                }
                else{
                    print $emblout $saved_line;
                    $printNext=1;
                    $saved_line="";
                    $skip_line=1;
                    print $line;
                }
            }
        }
    }

    if($printNext and !$skip_line){
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

 		if($printNext  and !$skip_line){
    		print $emblout $line;
    	}
    }

}

print "On $nbIDTotal headers there are $nbIDfound that has been modified properly.\n";
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


=item B<--bac>

Bolean. Believe in AC line. Instead of looking at the sequence size, look at the AC line (the second one of each record) as common information for the two files.

=item  B<--out>, B<--output> or B<-o>

The output will be the EMBL file with the record "headers" modified

=item B<--help> or B<-h>

Display this helpful text.

=back

=cut
