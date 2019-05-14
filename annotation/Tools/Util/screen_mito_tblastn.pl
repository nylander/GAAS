#!/usr/bin/env perl

use strict;
use warnings;
use IO::File;
use Bio::DB::Fasta;
use File::Basename;
use Getopt::Long;
use Pod::Usage;

my @copyARGV=@ARGV;

my $opt_output = undef;
my $tabfile = undef;
my $opt_genome= undef;
my $help= undef;

my $header = qq{
########################################################
# NBIS 2018 - Sweden                                   #  
# jacques.dainat\@nbis.se                              #
# Please cite NBIS (www.nbis.se) when using this tool. #
########################################################
};

if ( !GetOptions("tab=s"		  => \$tabfile,
		    	"o|out=s" 	  => \$opt_output,
                "g|genome=s"  => \$opt_genome,
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

### output
my $ostream     = IO::File->new();

# Manage Output
if(defined($opt_output))
{
$ostream->open( $opt_output, 'w' ) or
  croak(
     sprintf( "Can not open '%s' for reading: %s", $opt_output, $! ) );
}
else{
  $ostream->fdopen( fileno(STDOUT), 'w' ) or
      croak( sprintf( "Can not open STDOUT for writing: %s", $! ) );
}



#####
my %allIDs; # save ID in lower case to avoid cast problems
my $db = undef;
if ($opt_genome){
    my $nbFastaSeq=0;
    $db = Bio::DB::Fasta->new($opt_genome);
    my @ids      = $db->get_all_primary_ids;
    foreach my $id (@ids ){$allIDs{lc($id)}=$id;}
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
        my $ID = $list[1];
        my $mito_gene = $list[0];
        $mito_gene  =~ /\w+\|\w+\|([^_]+).*/;
        $mito_gene = $1;
        #print " my mytogene = $mito_gene \n";
        my $start=undef;
        my $end= undef;
        if($list[8]< $list[9]){
            $start = $list[8];
            $end = $list[9];
        }
        else{
            $start = $list[9];
            $end = $list[8];
        }       
         #print $start." ".$end."\n";
        push (@{$info{$ID}}, [$start, $end, $mito_gene]);
}
    
my %omni;
my %nbMitoGeneByContig;

foreach my $contig (keys %info){

    my @uniq_list;
    #print "contig: ".$contig."\n";

    my $prev_start = -1;
    my $prev_end = -1;
    my $start = -1;
    my $end = -1;
    my $printed=undef;

    foreach my $tuple (sort {$a->[0] <=> $b->[0] }  @{$info{$contig}}){

        $nbMitoGeneByContig{$contig}{$tuple->[2]}++;
        #print $tuple->[2];exit;
       
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

    if ($opt_genome){
        print $ostream "SequenceID\tNumber_of_non_ovelaping_Hit\tNb_mito_gene\tTotal_hit_size\tSize_sequence\t%_Sequence_covered_by_hit\tGene_names\n";
        # sort by number of non-overlaping hits
        foreach my $contig (sort { @{$omni{$a}} <=> @{$omni{$b}} } keys %omni){
            
            #compute length of the contig
            my $seq_id_correct = $allIDs{lc($contig)};
            my $seq     = $db->get_Seq_by_id($seq_id_correct);
            my $length  = $seq->length;
            
            #compute % of seq covered by mito hits
            my $goodxGenome=sprintf("%0.2f",($size{$contig}*100)/$length);
            
            my $nbMitoGene = keys %{$nbMitoGeneByContig{$contig}};
            
            my @geneList=();
            foreach my $key (sort keys %{$nbMitoGeneByContig{$contig}}){
             push @geneList, $key   
            }
            print $ostream $contig."\t".@{$omni{$contig}}."\t".$nbMitoGene."\t".$size{$contig}."\t".$length."\t".$goodxGenome."\t".join(",", @geneList)."\n";
        }

    }
    else{
        print $ostream "SequenceID\tNumber_of_non_ovelaping_Hit\tNb_mito_gene\tTotal_hit_size\tGene_names\n";
        # sort by number of non-overlaping hits
        foreach my $contig (sort { @{$omni{$a}} <=> @{$omni{$b}} } keys %omni){
            my $nbMitoGene = keys %{$nbMitoGeneByContig{$contig}};
            
            my @geneList=();
            foreach my $key (sort keys %{$nbMitoGeneByContig{$contig}}){
             push @geneList, $key
            }
            
            print $ostream $contig."\t".@{$omni{$contig}}."\t".$nbMitoGene."\t".$size{$contig}."\t".join(",", @geneList)."\n";
        }
    }

__END__


=head1 NAME
 
Based on a default blast tabulated  output ( -outfmt 6 => qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore) the script aims to tell you 
for each sequence of your assembly, how many non-overlaping mito hits have been found, the number of mito genes that have a hit and the total size in bp of those hits 
(overlaping part counted only once). When the assembly is provided, 2 new columns are displayed, the sie of the Sequence and the % part covered by mito hits. 
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

Input tabulated blast file -outfmt 6

=item  B<--out>, B<--output> or B<-o>

The output will be the EMBL file with the record "headers" modified

=item  B<--genome> or B<-g>

Optional. Genome in fasta format. Allow to calculate the mapping coverage.

=item B<--help> or B<-h>

Display this helpful text.

=back

=cut
