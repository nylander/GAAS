#!/usr/bin/env perl -w

## BILS 2015
## jacques.dainat@bils.se

## TO DO => Deal With sequences. Write the DNA sequence of the "source" primary tag within the output gff3  


use strict;
use Pod::Usage;
use Getopt::Long;
use Bio::Tools::GFF;
use Data::Dumper;
use Bio::SeqIO;

my $header = qq{
########################################################
# BILS 2015 - Sweden                                   #  
# jacques.dainat\@bils.se                               #
# Please cite BILS (www.bils.se) when using this tool. #
########################################################
};

my $outfile = undef;
my $gb = undef;
my $primaryTags = undef;
my $discard = undef;
my $keep = undef;
my $help;

if( !GetOptions(
    "help" => \$help,
    "gb=s" => \$gb,
    "ptag|t=s" => \$primaryTags,
    "d|s" => \$discard,
    "k" => \$keep,
    "outfile|output|o|out|gff=s" => \$outfile))
{
    pod2usage( { -message => "Failed to parse command line\n$header",
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($help) {
    pod2usage( { -verbose => 2,
                 -exitval => 2,
                 -message => "$header\n" } );
}

if ( ! (defined($gb)) ){
    pod2usage( {
           -message => "$header\nMissing the --gb argument",
           -verbose => 0,
           -exitval => 1 } );
}

##################
# MANAGE OPTION  #
if($discard and $keep){
  print "Cannot discard and keep the same primary tag. You have to choose if you want to discard it or to keep it.\n";
}

### If primaryTags given, parse them:

my @listprimaryTags;
if ($primaryTags){
  @listprimaryTags= split(/,/, $primaryTags);

  if($discard){ 
    print "We will not keep the following primary tag:\n";
    foreach my $tag (@listprimaryTags){ 
      print $tag,"\n";
    }
  }
  elsif($keep){ # Attribute we have to replace by a new name
    print "We will keep only the following primary tag:\n";
    foreach my $tag (@listprimaryTags){ 
      print $tag,"\n";
    }
  }
  else{print "You gave a list of primary tag wihtout telling me what you want I do with. Discard them or keep only them ?\n";}
}


##################
# MANAGE OUTPUT  #
my $gff_out;
if ($outfile) {
open(my $fh, '>', $outfile) or die "Could not open file '$outfile' $!";
  $gff_out= Bio::Tools::GFF->new(-fh => $fh, -gff_version => 3);
}
else{
  $gff_out = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 3);
}

### Read gb input file. 
my $gb_in = Bio::SeqIO->new(-file => $gb, -format => 'genbank');


### MAIN ###

while( my $seq_obj = $gb_in->next_seq) {

  for my $feat_obj ($seq_obj->get_SeqFeatures) {          
    my $skipit=undef;

    # In case we should discard some
    if($discard){
      
      foreach my $pTag (@listprimaryTags){
        if(lc($pTag) eq lc($feat_obj->primary_tag)){
          $skipit=1;last;
        }
      }
    }
    # In case we should keep only some
    elsif($keep){
      my $skipit=1;
      foreach my $pTag (@listprimaryTags){
        if(lc($pTag) eq lc($feat_obj->primary_tag)){
          $skipit=undef;last;
        }
      }

    }
    
    if(! $skipit){
        $gff_out->write_feature($feat_obj);
    }   
}
  #my $seqObject = Bio::Seq->new(-seq => $db->seq($seq_id));

  #$gff_out->write_feature($seq);
}

__END__

=head1 NAME

gb2gff.pl -
The script take a GeneBank file as input, and will translate it in gff format.

=head1 SYNOPSIS

    ./gb2gff.pl --gb infile.gb [ -o outfile ]

=head1 OPTIONS

=over 8

=item B<--gb>

Input EMBL file that will be read 

=item B<-primary_tag>, B<--pt>, B<-t>

List of "primary tag". Useful to discard or keep specific features.
The tags have to be separated by a coma.

=item B<-d>

Means that primary tags provided by the option "prinary_tag" will be discarded.

=item B<-d>

Means that only primary tags provided by the option "primary_tag" will be kept.

=item B<-o> , B<--output> , B<--out> , B<--outfile> or B<--gff>

Output GFF file. If no output file is specified, the output will be
written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.

=back

=cut