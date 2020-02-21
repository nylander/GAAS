#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use List::MoreUtils qw(uniq);
use Pod::Usage;
use Bio::Tools::GFF;
use IO::File;
use GAAS::GFF3::Omniscient;

my $header = qq{
########################################################
# NBIS 2015 - Sweden                                   #
# jacques.dainat\@nbis.se                               #
# Please cite NBIS (www.nbis.se) when using this tool. #
########################################################
};

my $opt_output= undef;
my $opt_score = undef;
my $opt_test = undef;
my $opt_gff = undef;
my $opt_help;

# OPTION MANAGMENT
my @copyARGV=@ARGV;
if ( !GetOptions( 'f|ref|reffile|gff=s' => \$opt_gff,
                  's|score|v=f'         => \$opt_score,
                  't|test=s'            => \$opt_test,
                  'o|output=s'          => \$opt_output,
                  'h|help!'             => \$opt_help ) )
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

if ( ! $opt_gff or ! defined($opt_score) or ! $opt_test ){
    pod2usage( {
           -message => "$header\nAt least 3 parameters are mandatory:\n1) Input reference gff file (--f)\n".
           "2) score use for filtering (between 0 and 1) with option --v\n3) test to apply (> < = >= <=) with the option -t. And don't forget to quote you parameter if it contains the character < or >.\n\n",
           -verbose => 0,
           -exitval => 1 } );
}

###############
# Test options
if($opt_score < 0 or $opt_score > 1){
  print "The value of the score option is Wrong: $opt_score.\n We want a value between 0 and 1.";exit;
}
if($opt_test ne "<" and $opt_test ne ">" and $opt_test ne "<=" and $opt_test ne ">=" and $opt_test ne "="){
  print "The test to apply is Wrong: $opt_test.\nWe want something among this list: <,>,<=,>= or =.";exit;
}

###############
# Manage Output

# FOR REPORT
my $ostreamReport = \*STDOUT or die ( sprintf( "Can not open '%s' for writing %s", "STDOUT", $! ));

## FOR GFF FILE
my $ostream_ok  = IO::File->new();
my $ostream_discarded  = IO::File->new();
if($opt_output){
  print $opt_output;
  $opt_output =~ s/\.gff$//g;
  print "after ".$opt_output;
  my $out_ok = "$opt_output".".gff";
  my $out_discarded = $opt_output."_discarded.gff";

  if(-f $out_ok or -f $out_discarded){
    print "File already exist.\n";exit;
  }
  $ostream_ok->open( $out_ok, 'w' ) or croak( sprintf( "Can not open '%s' for reading: %s", $opt_output, $! ) );
  $ostream_discarded->open( $out_discarded, 'w' ) or croak( sprintf( "Can not open '%s' for reading: %s", $opt_output, $! ) );
}
else{
  $ostream_ok->fdopen( fileno(STDOUT), 'w' ) or
      croak( sprintf( "Can not open STDOUT for writing: %s", $! ) );
  $ostream_discarded->fdopen( fileno(STDOUT), 'w' ) or
      croak( sprintf( "Can not open STDOUT for writing: %s", $! ) );
}
my $gffout_ok = Bio::Tools::GFF->new( -fh => $ostream_ok , -gff_version => 3) or croak( sprintf( "Can not open STDOUT for writing: %s", $! ) );
my $gffout_discarded = Bio::Tools::GFF->new( -fh => $ostream_discarded , -gff_version => 3) or croak( sprintf( "Can not open STDOUT for writing: %s", $! ) );



                                                      #######################
                                                      #        MAIN         #
#                     >>>>>>>>>>>>>>>>>>>>>>>>>       #######################       <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# start with some interesting information
my $stringPrint = strftime "%m/%d/%Y at %Hh%Mm%Ss", localtime;
$stringPrint .= "\nusage: $0 @copyARGV\n";

######################
### Parse GFF input #
my ($hash_omniscient, $hash_mRNAGeneLink) =  slurp_gff3_file_JD({
                                                               input => $opt_gff,
                                                               verbose => 1
                                                               });
print("Parsing Finished\n\n");
### END Parse GFF input #
#########################


###########################
# Main compute
my @listIDl2discarded;
my @listIDl2ok;
foreach my $tag_l1 (keys %{$hash_omniscient->{'level1'}}){
  foreach my $id_l1 (keys %{$hash_omniscient->{'level1'}{$tag_l1}}){

    #################
    # == LEVEL 2 == #
    #################
    foreach my $tag_l2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
      if ( exists ($hash_omniscient->{'level2'}{$tag_l2}{$id_l1} ) ){
        foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}}) {
          if($feature_level2->has_tag('_AED')){
            my $AED_score=$feature_level2->_tag_value('_AED');
            my $id_level2=lc($feature_level2->_tag_value('ID'));

            if ($opt_test eq ">"){
              if ($AED_score > $opt_score){
                push @listIDl2ok, $id_level2 ;
              }else{push @listIDl2discarded, $id_level2 ;}
            }
            if ($opt_test eq "<"){
              if ($AED_score < $opt_score){
                push @listIDl2ok, $id_level2 ;
              }else{push @listIDl2discarded, $id_level2 ;}
            }
            if ($opt_test eq "="){
              if ($AED_score == $opt_score){
                push @listIDl2ok, $id_level2 ;
              }else{push @listIDl2discarded, $id_level2 ;}
            }
            if ($opt_test eq "<="){
              if ($AED_score <= $opt_score){
                push @listIDl2ok, $id_level2 ;
              }else{push @listIDl2discarded, $id_level2 ;}
            }
            if ($opt_test eq ">="){
              if ($AED_score >= $opt_score){
                push @listIDl2ok, $id_level2 ;
              }else{push @listIDl2discarded, $id_level2 ;}
            }
          }
          else{
            print "WARNING: _AED attribute not found for feature ".$ostreamReport->write_feature($feature_level2);
          }
        }
      }
    }
  }
}

# remove duplicate in case several option tends to give the same case
if(@listIDl2ok){
  my $sizeList= @listIDl2ok;
  $stringPrint.= "$sizeList RNA(s) that reach your quality request. ($opt_test $opt_score)\n";
  my @listIDl2okUniq = uniq(@listIDl2ok);
  my $omniscient_ok = create_omniscient_from_idlevel2list($hash_omniscient, $hash_mRNAGeneLink, \@listIDl2okUniq);
  print_omniscient($omniscient_ok, $gffout_ok);
}
if(@listIDl2discarded){
  my $sizeList= @listIDl2discarded;
  $stringPrint.= "$sizeList RNA(s) discarded because don't reach your quality request. ($opt_test $opt_score)\n";
  my @listIDl2discardedUniq = uniq(@listIDl2discarded);
  my $omniscient_discarded = create_omniscient_from_idlevel2list($hash_omniscient, $hash_mRNAGeneLink, \@listIDl2discarded);
  print_omniscient($omniscient_discarded, $gffout_discarded);
}

#Print Info Output
print $ostreamReport $stringPrint;


=head1 NAME

maker_select_models_by_AED_score.pl -
The script take a gff3 file as input. -

The result is written to the specified output file, or to STDOUT.
Remark: If there is duplicate in the file they will be removed in the output. In that case you should be informed.

=head1 SYNOPSIS

    ./maker_select_models_by_AED_score.pl -f infile.gff -v 1 -t = [ --output outfile ]
    ./maker_select_models_by_AED_score.pl --help

=head1 OPTIONS

=over 8

=item B<-f>, B<--reffile>, B<--gff>  or B<-ref>

Input GFF3 file that will be read

=item B<-v>, B<--score> or B<-s>

Score use for filtering (between 0 and 1) with option.
Can be a float.

=item B<-t> or B<--test>
Test to apply (> < = >= <=). If you us one of these two character >, <, please don't forget to quote you parameter liket that "<=". Else your terminal will complain.

=item B<-o> or B<--output>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.

=back

=cut
