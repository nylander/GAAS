#!/usr/bin/perl

use strict;
use Env qw (@B2G4PIPEPATH);
use Getopt::Long;
use Scalar::Util qw(openhandle);
use Time::Piece;
use Time::Seconds;

my $usage = qq{
perl run_blast2go.pl
  Getting help:
    [--help]

  Input:
    [--blast filename]
		The name of the blast file  (-m 7)
    [--ipr file]
		The name of the iprscan xml file
  Ouput:    
    [--output filename]
        The name of the output file. 

  Misc:
	[--mem amount]
		Specify the memory to be used by Blast2Go (e.g. 800m, 10G)
		Default: 1G

};

my $blast = undef;
my $ipr = undef;
my $output = undef;
my $mem = "1G";

my $help;

GetOptions(
    "help" => \$help,
    "blast=s" => \$blast,
    "ipr=s" => \$ipr,
	"mem=s" => \$mem,
    "output=s" => \$output);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}


if (-f "@B2G4PIPEPATH/b2gPipe.properties") {
	print "Found Blast2Go installation, starting analysis\n";
} else {
	die "Couldn't find the Blast2Go binary - make sure you load the relevant module!\n"
}

# $B2G4PIPE/*:$B2G4PIPE/ext/*:

runcmd("java -Xmx$mem -cp @B2G4PIPEPATH/*:@B2G4PIPEPATH/ext/* es.blast2go.prog.B2GAnnotPipe -in $blast -out $output -prop @B2G4PIPEPATH/b2gPipe.properties -ips $ipr -annot");


sub msg {
  my $t = localtime;
  my $line = "[".$t->hms."] @_\n";
}

sub runcmd {
  msg("Running:", @_);
  system(@_)==0 or err("Could not run command:", @_);
}

sub err {
  msg(@_);
  exit(2);
}


