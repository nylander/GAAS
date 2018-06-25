#!/usr/bin/env perl

use strict;
use warnings;
use IO::File;
use File::Basename;
use Getopt::Long;
use Pod::Usage;

my @copyARGV=@ARGV;

my $opt_output = undef;
my $tabfile = undef;
my $help= undef;

my $header = qq{
########################################################
# BILS 2016 - Sweden                                   #  
# jacques.dainat\@bils.se                               #
# Please cite BILS (www.bils.se) when using this tool. #
########################################################
};

if ( !GetOptions("tab=s"		  => \$tabfile,
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

if ( ! defined( $tabfile) ) {
    pod2usage( {
           -message => "$header\nMust specify at least 1 file. \n",
           -verbose => 0,
           -exitval => 1 } );
}

##### Stream in 1 
my $fh1;
if ($tabfile) {
  open($fh1, '<', $tabfile) or die "Could not open file '$tabfile' $!";
}


my %info;

while( my $line = <$fh1>)  {   
    
    if( $line =~ m/^#/){ next; }

    my @list = split(/\s/,$line);
    my $ID = $list[0];
    my $start = $list[9];
    my $end = $list[10];
	#print $start." ".$end."\n";
    push (@{$info{$ID}}, [$start, $end]);

}
    
my %omni;
foreach my $contig (keys %info){

    my @uniq_list;
    #print "contig: ".$contig."\n";

    my $prev_start = -1;
    my $prev_end = -1;
    my $start = -1;
    my $end = -1;
    my $printed=undef;

    foreach my $tuple (sort {$a->[0] <=> $b->[0] }  @{$info{$contig}}){
       
        $start = @$tuple[0];
        $end = @$tuple[1];

        if ( ($prev_start <= $end) and ($prev_end >= $start) ){ #it overlaps or are consecutive
            #print "it overlaps\n";
            if ($end > $prev_end){
                $prev_end = $end;
            }  
        } 
        elsif($start > $prev_end){
            if($prev_start != -1){
                push (@uniq_list, [$prev_start,$prev_end]);
                #print "I push the tuple [$prev_start,$prev_end]\n";
            }
            $prev_start = $start ;
            $prev_end = $end ;                
        }

    }
    # Deal with the last round
    push (@uniq_list, [$prev_start,$prev_end]);
    #print "I push the last tuple [$prev_start,$prev_end]\n";

    push (@{$omni{$contig}}, @uniq_list)
}
    
    #calculate bp incremented non-overlaping hit size
    my %size;
    foreach my $contig (keys %omni){
        foreach my $tuple ( @{$omni{$contig}} ){
            $size{$contig}+=(@$tuple[1]-@$tuple[0]+1);
        }
    }

    print "SequenceID\tNumber_of_Hit\tTotal_bp\n";
    # sort by number of non-overlaping hits
    foreach my $contig (sort { @{$omni{$a}} <=> @{$omni{$b}} } keys %omni){
        print $contig."\t".@{$omni{$contig}}."\t".$size{$contig}."\n";

    }


__END__


=head1 NAME
 
Based on a balst output ( -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle') the script aims to tell you how many non-overlaping hit has been found
by sequence. One hit may roughly be considered as one gene. It gives also the total size of those hits in bp.
The script aims to help determining the contigs from an assembly which are mitochondrial. An assembly graph could be helpful to check if the suspicious (those that might be mitochondrial) contigs sounds to be circular
as expected for a mitochondrial genome.

=head1 SYNOPSIS

    ./screen_mito.pl --tab=infile -o=outFile 
    ./screen_mito.pl --help

Mitochondrial genome size (from wikipedia)

Genome Type Kingdom Introns Size    Shape   Description
1   Animal  No  11–28kbp    Circular    Single molecule
2   Fungi, Plant, Protista  Yes 19–1000kbp  Circular    Single molecule
3   Fungi, Plant, Protista  No  20–1000kbp  Circular    Large molecule and small plasmid like structures
4   Protista    No  1–200kbp    Circular    Heterogeneous group of molecules
5   Fungi, Plant, Protista  No  1–200kbp    Linear  Homogeneous group of molecules
6   Protista    No  1–200kbp    Linear  Heterogeneous group of molecules



=head1 OPTIONS

=over 8

=item B<--tab>

Input tabulated blast file -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle'

=item  B<--out>, B<--output> or B<-o>

The output will be the EMBL file with the record "headers" modified

=item B<--help> or B<-h>

Display this helpful text.

=back

=cut
