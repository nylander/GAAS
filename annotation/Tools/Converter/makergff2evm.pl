#!/usr/bin/env perl

# A filter for Uniprot and RefSeq fasta files that makes the fasta
# headers a bit more terse.  Reads from STDIN, writes to STDOUT.
#

use strict;
use warnings;

my $first_line = <STDIN>;
chomp($first_line);


my %parsers = (
    'null' => sub {
        my ($line) = @_;
        return $line;
    },

    'est' => sub {
        my ($line) = @_;
			
		$line =~ s/match_part/EST_match/g ;
			
		return $line;

    },
    'protein' => sub {
        my ($line) = @_;

        $line =~ s/match_part/nucleotide_to_protein_match/g ;

        return $line;
    } );

my $parser = 'null';

if    ( $first_line =~ /^.*est2genome.*/ ) { $parser = 'est'; }
elsif ( $first_line =~ /^.*Cufflinks.*/ ) { $parser = 'est'; }
elsif ( $first_line =~ /^.*protein2genome.*/ ) { $parser = 'protein'; }

print( $parsers{$parser}($first_line), "\n" );

while ( my $line = <STDIN> ) {
    chomp($line);
    next unless ($line =~ /match_part/ );
    print $parsers{$parser}($line), "\n";
}
