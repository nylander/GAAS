#!/usr/bin/env perl

use Carp;
use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use IO::File ;
use Bio::Tools::GFF;
use BILS::Handler::GFF3handler qw(:Ok);
use BILS::Handler::GXFhandler qw(:Ok);

my $start_run = time();
my $header = qq{
########################################################
# BILS 2018 - Sweden                                   #
# jacques.dainat\@nbis.se                               #
# Please cite NBIS (www.nbis.se) when using this tool. #
########################################################
};

my $inputFile=undef;
my $outfile=undef;
my $outformat=undef;
my $primaryTag=undef;
my $opt_help = 0;
my $locus_tag=undef;
my $tag_in="ID";

Getopt::Long::Configure ('bundling');
if ( !GetOptions ('file|input|gff=s' => \$inputFile,
      'to|lo=s' => \$locus_tag,
      'ti|li=s' => \$tag_in,
      "p|type|l=s" => \$primaryTag,
      'o|output=s' => \$outfile,
      'of=i' => \$outformat,
      'h|help!'         => \$opt_help )  )
{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

if ($opt_help) {
    pod2usage( { -verbose => 2,
                 -exitval => 2,
                 -message => "$header\n" } );
}

if ((!defined($inputFile)) ){
   pod2usage( { -message => "$header\nAt least 1 parameter is mandatory: -i",
                 -verbose => 0,
                 -exitval => 1 } );
}

my $ostream     = IO::File->new();

# Manage input fasta file
my $format = select_gff_format($inputFile);
my $ref_in = Bio::Tools::GFF->new(-file => $inputFile, -gff_version => $format);

# Manage Output
if(! $outformat){
  $outformat=$format;
}

# Manage Output
my $gffout;
if ($outfile) {
  $outfile=~ s/.gff//g;
  open(my $fh, '>', $outfile.".gff") or die "Could not open file '$outfile' $!";
  $gffout= Bio::Tools::GFF->new(-fh => $fh, -gff_version => $outformat );

}
else{
  $gffout = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => $outformat );
}

#define the locus tag
if(! $locus_tag){
  $locus_tag="locus_tag";
}

# Manage $primaryTag
my @ptagList;
if(! $primaryTag or $primaryTag eq "all"){
  print "We will work on attributes from all features\n";
  push(@ptagList, "all");
}elsif($primaryTag =~/^level[123]$/){
  print "We will work on attributes from all the $primaryTag features\n";
  push(@ptagList, $primaryTag);
}else{
   @ptagList= split(/,/, $primaryTag);
   foreach my $tag (@ptagList){
      if($tag =~/^level[123]$/){
        print "We will work on attributes from all the $tag features\n";
      }
      else{
       print "We will work on attributes from $tag feature.\n";
      }
   }
}

#time to calcul progression
my $startP=time;
my $nbLine=`wc -l < $inputFile`;
$nbLine =~ s/ //g;
chomp $nbLine;
print "$nbLine line to process...\n";

my $geneName=undef;
my $line_cpt=0;
while (my $feature = $ref_in->next_feature() ) {
  $line_cpt++;


    manage_attributes($feature, \@ptagList, $locus_tag);
    $gffout->write_feature($feature);

  #Display progression
  if ((30 - (time - $startP)) < 0) {
    my $done = ($line_cpt*100)/$nbLine;
    $done = sprintf ('%.0f', $done);
        print "\rProgression : $done % processed.\n";
    $startP= time;
  }
}


##Last round
my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job done in $run_time seconds\n";


#######################################################################################################################
        ####################
         #     methods    #
          ################
           ##############
            ############
             ##########
              ########
               ######
                ####
                 ##

sub  manage_attributes{
  my  ($feature, $ptagList, $locus_tag)=@_;

  if($feature->has_tag($tag_in)){
    print "Ihave\n";
    my $level = get_level($feature);
    my $primary_tag=$feature->primary_tag;
    my $id=$feature->_tag_value($tag_in);

    # check primary tag (feature type) to handle
    foreach my $ptag (@$ptagList){
      if($ptag eq "all"){
        print "that $locus_tag $id\n";
        create_or_replace_tag($feature,$locus_tag, $id);
      }
      elsif(lc($ptag) eq $level){
        create_or_replace_tag($feature,$locus_tag, $id);
      }
      elsif(lc($ptag) eq lc($primary_tag) ){
        create_or_replace_tag($feature,$locus_tag, $id);
      }
    }
  }
}

__END__

=head1 NAME

gff3_add_locus_tag.pl -
add a locus tag based on the ID of the a feature (by defaul gene feature).

=head1 SYNOPSIS

    gff3_add_locus_tag.pl --gff <input file> [-o <output file>]
    gff3_add_locus_tag.pl --help

=head1 OPTIONS

=over 8

=item B<--gff>, B<--file> or B<--input>

STRING: Input gff file that will be read.

=item B<-p>,  B<--type> or  B<-l>

primary tag option, case insensitive, list. Allow to specied the feature types that will be handled.
You can specified a specific feature by given its primary tag name (column 3) as: cds, Gene, MrNa
You can specify directly all the feature of a particular level:
      level2=mRNA,ncRNA,tRNA,etc
      level3=CDS,exon,UTR,etc
By default all feature are taking in account. fill the option by the value "all" will have the same behaviour.

=item B<--lo> or B<--tp>

Locus tag output, by defaut it will be called locus_tag, but using this option you can specied the name of this attribute.

=item B<--li> or B<--ti>

Tag input, by default the value of the attribute ID is used to fill the locus tag output, but you can chose whichever you want with this option.

=item B<--of>

Output format, if no ouput format is given, the same as the input one detected will be used. Otherwise you can force to have a gff version 1 or 2 or 3 by giving the corresponding number.

=item B<-o> or B<--output>

STRING: Output file.  If no output file is specified, the output will be written to STDOUT. The result is in tabulate format.

=item B<--help> or B<-h>

Display this helpful text.

=back

=cut
