#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case bundling);
use IO::File;
use File::Basename;

my $EXPECTED_NUM_COLUMNS = 9;
my $in = \*STDIN;
my $out_dir;

parse_options();
split_gff3();

sub parse_options {
	my $help;
	my $infile;
	GetOptions("i|input=s"		  => \$infile,
		   "d|output_directory=s" => \$out_dir,
		   "h|help"		  => \$help);
	print_usage() if $help;
	die "Missing required output directory (-d) option\n" if !$out_dir;
	$in = new IO::File($infile) or die "Error reading gff3: $!\n"
		if $infile;
}

sub print_usage {
	my $progname = basename($0);
	die << "END";
usage: $progname
	[-i|--input <gff3_file>]
	-d|--output_directory <directory_to_output_to>
	[-h|--help]

	i: input GFF3 [default: stdin]
	d: directory to write split fasta files into
END
}

sub split_gff3 {
	my %data = ();

	while (my $line = <$in>) {
		chomp $line;
		my @tokens = split /\t+/, $line;
		next if scalar(@tokens) != $EXPECTED_NUM_COLUMNS;
		next if $tokens[1] eq ".";
		push @{$data{$tokens[1]}}, \@tokens;
	}

	my @sources = keys(%data);

	for (my $i = 0; $i < scalar(@sources); ++$i) {
		my $source = $sources[$i];
		my $out = new IO::File(">$out_dir/$source.gff") or
			die "Error writing $out_dir/$source.gff: $!\n";
		print $out "##gff-version 3\n";
		foreach my $gff (@{$data{$source}}) {
			print $out join("\t", @{$gff}), "\n";
		}
		$out->close();
	}
}

