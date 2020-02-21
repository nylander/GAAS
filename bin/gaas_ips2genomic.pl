#!/usr/bin/env perl

#-----------------------------------------------------------------------
# Reads mRNA and CDS data from a GFF file and then maps the output from
# InterProScan to genomic locations given InterProScan results in TSV
# format.
#

use strict;
use warnings;
use Getopt::Long;
use IO::File;

my $usage = qq{
perl my_script.pl
  Getting help:
    [--help]
  List of all available databases
	[--gff_file] 
	A GFF file with gene models (reference data)
	[--ips_file]
	An Interpro TSV formatted annotation
};

my $gff_file = undef;
my $ips_file = undef;
my $help;

GetOptions(
	"help" => \$help,
	"gff_file=s" => \$gff_file,
	"ips_file=s" => \$ips_file);	

sub parse_gff_line
{
    my ($line) = @_;

    my @fields = split( /\t/, $line );

    my %entry;

    $entry{'seqid'}  = $fields[0];
    $entry{'type'}   = $fields[2];
    $entry{'start'}  = $fields[3];
    $entry{'end'}    = $fields[4];
    $entry{'strand'} = $fields[6];

    $fields[8] =~ /ID=([^;]+)/;
    $entry{'id'} = $1;

    $fields[8] =~ /Parent=([^;]+)/;
    $entry{'parent'} = $1;

    return \%entry;
}

my $gff_in = IO::File->new($gff_file)
  or die( sprintf( "Unable to open GFF file '%s' for reading:\n%s\n",
                   $gff_file, $! ) );

my $ips_in = IO::File->new($ips_file)
  or die( sprintf( "Unable to open "
                     . "InterProScan result file '%s' for reading:\n%s\n",
                   $ips_file, $! ) );

my %mRNA;

print( STDERR ">> Reading GFF file...\n" );
while ( my $line = $gff_in->getline() ) {
    chomp($line);

    my $entry = parse_gff_line($line);

    # Only care about the 'mRNA' and 'CDS' entries.  Hook all CDS
    # entries up to the corresponding mRNA.
    if ( $entry->{'type'} eq 'mRNA' ) {
        $mRNA{ $entry->{'id'} } = $entry;
    }
    elsif ( $entry->{'type'} eq 'CDS' ) {
        if ( exists( $mRNA{ $entry->{'parent'} } ) ) {
            push( @{ $mRNA{ $entry->{'parent'} }{'CDS'} }, $entry );
        }
        else {
            die( sprintf( "No parent '%s' for CDS '%s'",
                          $entry->{'parent'}, $entry->{'id'} ) );
        }
    }
}

$gff_in->close();

print( STDERR ">> Calculating CDS coordinates...\n" );
foreach my $mRNA_entry ( values(%mRNA) ) {
    my $cds_length_sum = 0;

    $mRNA_entry->{'CDS'} = [ sort { $a->{'start'} <=> $b->{'start'} }
                             @{ $mRNA_entry->{'CDS'} } ];

    foreach my $cds_entry ( @{ $mRNA_entry->{'CDS'} } ) {
        my $cds_length = $cds_entry->{'end'} - $cds_entry->{'start'} + 1;

        $cds_entry->{'CDS_start'} = $cds_length_sum + 1;
        $cds_entry->{'CDS_end'}   = $cds_length_sum + $cds_length;

        $cds_length_sum += $cds_length;
    }
}

print( STDERR ">> Processing InterProScan results...\n" );
while ( my $line = $ips_in->getline() ) {
    chomp($line);

    my @ips_fields = split( /\t/, $line );

    my $mRNA_id = $ips_fields[0];
    if ( !exists( $mRNA{$mRNA_id} ) ) {
        warn(
            sprintf( "mRNA ID '%s' not found in GFF file, skipping\n", $mRNA_id ) );
        next;
    }

    my $mRNA_entry = $mRNA{$mRNA_id};

    my $hit_cds_start = 3 * ( $ips_fields[6] - 1 ) + 1;
    my $hit_cds_end   = 3 * ( $ips_fields[7] - 1 ) + 1;

    my $hit_start;
    my $hit_end;
    my @hit_coords;    # List of "match_part" coordinates.

    # Go through the CDS entries for this mRNA until we've found the
    # correct one for both hit start and hit end.
    foreach my $cds_entry ( @{ $mRNA_entry->{'CDS'} } ) {
        if (    $hit_cds_start >= $cds_entry->{'CDS_start'}
             && $hit_cds_start <= $cds_entry->{'CDS_end'} )
        {
            # Start of hit is in this CDS.
            $hit_start = $cds_entry->{'start'} +
              ( $hit_cds_start - $cds_entry->{'CDS_start'} );

            # This "match_part" starts part way into this CDS.  Its end
            # is still unknown.
            @hit_coords = ( [ $hit_start, undef ] );
        }
        elsif ( defined($hit_start) ) {

            # CDS for hit start has been found already and this CDS is
            # part of the hit, so this "match_part" starts at the start
            # of the CDS.  Its end is still unknown.
            push( @hit_coords, [ $cds_entry->{'start'}, undef ] );
        }

        if (    $hit_cds_end >= $cds_entry->{'CDS_start'}
             && $hit_cds_end <= $cds_entry->{'CDS_end'} )
        {
            # End of hit is in this CDS.
            $hit_end =
              $cds_entry->{'start'} + ( $hit_cds_end - $cds_entry->{'CDS_start'} );

            # Complete the last "match_part" by filling in the hit end,
            # which is part way into this CDS.
            $hit_coords[-1][1] = $hit_end;

            last;    # Done with this protein match.
        }
        elsif ( defined($hit_start) ) {

            # This CDS is part of the hit, so this "match_part" ends at
            # the end of the CDS.
            $hit_coords[-1][1] = $cds_entry->{'end'};
        }
    } ## end foreach my $cds_entry ( @{ ...})

    if ( !defined($hit_start) || !defined($hit_end) ) {
        die( sprintf( "\nProtein match falls outside of CDS for %s::%s on %s\n",
                      $ips_fields[3], $ips_fields[4], $mRNA_id ) );
    }

    my $feature_id = sprintf( "%s:%s:%s",
                              $mRNA_entry->{'id'}, $ips_fields[3],
                              $ips_fields[4] );

    # Output the "protein_match" feature.
    printf( "%s\tinterproscan\tprotein_match\t%d\t%d\t%g\t%s\t.\t"
              . "ID=%s;"
              . "Name=%s:%s;"
              . "description=%s\n",
            $mRNA_entry->{'seqid'},  $hit_coords[0][0],
            $hit_coords[-1][1],      0,
            $mRNA_entry->{'strand'}, $feature_id,
            $ips_fields[3],          $ips_fields[4],
            $ips_fields[5] );

    # Output the "match_part" features.
    my $count = 0;
    foreach my $coords (@hit_coords) {
        printf( "%s\tinterproscan\tmatch_part\t%d\t%d\t%g\t%s\t.\t"
                  . "ID=%s:%d;"
                  . "Name=%s:%s;"
                  . "Parent=%s\n",
                $mRNA_entry->{'seqid'}, $coords->[0],            $coords->[1],
                $ips_fields[8],         $mRNA_entry->{'strand'}, $feature_id,
                ++$count,               $ips_fields[3],          $ips_fields[4],
                $feature_id );
    }

} ## end while ( my $line = $ips_in...)

$ips_in->close();

