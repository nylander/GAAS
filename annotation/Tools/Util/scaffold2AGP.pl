#!/usr/bin/env perl
# 
# Creates a AGP-file needed by e.g. EMBL for a scaffolded assembly
# 
# By Henrik Lantz, NBIS/Uppsala University, Sweden

# usage: perl scaffold2AGP.pl -i scaffoldfile.fasta -o scaffoldfile.agp

use warnings;
use strict;
use Bio::SeqIO ;
use Getopt::Long;

my $infile;
my $outfile;
my $contigcount=0;
my $usage = "\nUsage: perl scaffold2AGP.pl -i <scaffoldfile.fasta> -o <scaffoldfile.agp>\n";
$usage .= "Contigs will be saved in \"contigs.fasta\"\n";

GetOptions( 
            'i=s' => \$infile,     # scaffoldfile
            'o=s'   => \$outfile);

die "\nPlease provide filename(s)\n$usage\n" unless $infile;

open (AGP_FILE, ">$outfile") or
   die "\nPlease provide filename(s)\n$usage";


my $inseq = Bio::SeqIO->new('-file' => "<$infile",
               '-format' => 'Fasta' );
               
my $outseq = Bio::SeqIO->new(
            -file     => ">contigs.fasta",
            -format => 'fasta',
            );

#Read scaffolded FASTA-file
while (my $seq_obj = $inseq->next_seq ) {
  my $scaffold = $seq_obj->id;
  my $sequence = $seq_obj->seq;
  my $start=1;
  my $oldsum;
  my $newsum;
  my $count=0;
  my $rounded;
  
  next if ($scaffold =~ /^contig/i);
  foreach my $substring_sequence (split /(N{20,})/i, $sequence){
    my $type;
    my $substring_length = length($substring_sequence); 
    $count++;
    $oldsum=$start;  
    $newsum=$oldsum+$substring_length-1;
    
    if ($substring_sequence !~ m/^N+$/i){
      $type="W";  
      $contigcount++;
 	  $rounded=sprintf("%05s", $contigcount);
 	  my $contig_obj = Bio::Seq->new(-seq => "$substring_sequence",                        
                                           -display_id => "contig$rounded",                        
                                           -alphabet => "dna" );
      $outseq->write_seq($contig_obj);
    }  
    elsif ($substring_sequence =~ m/^N+$/i){
      $type="N";
    }
    $start += $substring_length;
    if ($type eq "W"){
     print AGP_FILE "$scaffold\t$oldsum\t$newsum\t$count\t$type\tcontig$rounded\t1\t$substring_length\t+\n";
    }
    if ($type eq "N"){
     print AGP_FILE "$scaffold\t$oldsum\t$newsum\t$count\t$type\t$substring_length\tscaffold\tyes\tpaired-ends\n";
    }
  }
}

close AGP_FILE;
