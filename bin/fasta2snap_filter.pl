#!/usr/bin/env perl
use strict;
use Bio::SeqIO;
 
my $filterfile = $ARGV[0];
my $seqs_to_keep = {};
open(my $FILTER, "<", $filterfile) or die("Error: could not open '$filterfile' $!");
while(<$FILTER>)
{
  chomp();
  $seqs_to_keep->{$_} = 1;
}
close($FILTER);
 
my $seqs = {};
my $loader = Bio::SeqIO->new(-fh => \*STDIN, -format => 'Fasta');
while(my $seq = $loader->next_seq)
{
  $seqs->{$seq->id} = $seq;
}
 
my @keys = sort(keys(%$seqs));
my $writer = Bio::SeqIO->new( -fh => \*STDOUT, -format => 'Fasta');
foreach my $seqid(@keys)
{
  $writer->write_seq($seqs->{$seqid}) if($seqs_to_keep->{$seqid});
}
