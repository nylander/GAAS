#!/usr/bin/env perl

# A filter for Uniprot and RefSeq fasta files that makes the fasta
# headers a bit more terse.  Reads from STDIN, writes to STDOUT.
#
## Note: Will pass any other fasta file unchanged.
## Note: For Uniprot, will also change any 'O' in the protein sequence
##       into 'K'.

use strict;
use warnings;

my $first_line = <STDIN>;
chomp($first_line);

if ( $first_line !~ /^>/ ) {
    die( sprintf( "This does not look like fasta formatted input:\n%s\n",
                  $first_line ) );
}

my %parsers = (
    'null' => sub {
        my ($line) = @_;
        return $line;
    },

    'uniprot' => sub {
        my ($line) = @_;


        if ( $line =~ /^>(?:sp|tr)\|([^|]+).*PE=(\d+) SV=(\d+)/ ) {
            return sprintf( ">%s.%d %d", $1, $3, $2 );
        }
        else {
            $line =~ tr/O/K/;
        }

        return $line;
    },
    'refseq' => sub {
        my ($line) = @_;

        if ( $line =~ /^>gi/ ) {
            return sprintf( ">%s", [ split( /\|/, $line ) ]->[3] );
        }

        return $line;
    } );

my $parser = 'null';

if    ( $first_line =~ /^>(?:sp|tr)/ ) { $parser = 'uniprot'; }
elsif ( $first_line =~ /^>gi/ ) { $parser = 'refseq'; }

print( $parsers{$parser}($first_line), "\n" );

while ( my $line = <STDIN> ) {
    chomp($line);
    print $parsers{$parser}($line), "\n";
}
