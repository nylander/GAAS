#!/usr/bin/perl

use strict;
use Pod::Usage;
use Getopt::Long;
use Bio::SeqIO ;
use Bio::DB::Fasta;
use Bio::Tools::GFF;
use BILS::Handler::GXFhandler qw(:Ok);
use BILS::Handler::GFF3handler qw(:Ok);


my $start_run = time();
my $codonTable=1;
my $opt_gfffile;
my $opt_fastafile;
my $opt_output;
my $opt_AA=undef;
my $opt_help = 0;
my $opt_extermityOnly=undef;
my $opt_upstreamRegion=undef;
my $opt_downRegion=undef;
my $opt_cdna=undef;
my $opt_type = 'cds';
my $width = 60; # line length printed

my $header = qq{
########################################################
# BILS 2016 - Sweden                                   #  
# jacques.dainat\@bils.se                               #
# Please cite BILS (www.bils.se) when using this tool. #
########################################################
};

# OPTION MANAGMENT
my @copyARGV=@ARGV;
if ( !GetOptions( 'g|gff=s' => \$opt_gfffile,
                  'f|fa|fasta=s' => \$opt_fastafile,
                  't=s' => \$opt_type,
                  'protein|p|aa' => \$opt_AA,
                  'cdna' => \$opt_cdna,
                  'ext|e' => \$opt_extermityOnly,
                  'table|codon|ct=i' => \$codonTable,
                  'up|5|five|upstream=i'      => \$opt_upstreamRegion,
                  'do|3|three|down|downstream=i'      => \$opt_downRegion,
                  'o|output=s'      => \$opt_output,
                  'h|help!'         => \$opt_help ) )
{
    pod2usage( { -message => "$header\nFailed to parse command line",
                 -verbose => 1,
                 -exitval => 1 } );
}
# shortcut for cdna
if($opt_cdna){$opt_type="exon";}

# Print Help and exit
if ($opt_help) {
    pod2usage( { -verbose => 2,
                 -exitval => 2,
                 -message => "$header \n" } );
}
 
if ( (! (defined($opt_gfffile)) ) or (! (defined($opt_fastafile)) ) ){
    pod2usage( {
           -message => "\nAt least 2 parametes are mandatory:\nInput reference gff file (-g);  Input reference fasta file (-f)\n\n".
           "Output is optional. Look at the help documentation to know more.\n",
           -verbose => 0,
           -exitval => 2 } );
}


my $ostream;
if ($opt_output) {
  $opt_output=~ s/.fasta//g;
  $opt_output=~ s/.fa//g;
  open(my $fh, '>', $opt_output.".fa") or die "Could not open file '$opt_output' $!";
  $ostream= Bio::SeqIO->new(-fh => $fh, -format => 'Fasta' );
}
else{
  $ostream = Bio::SeqIO->new(-fh => \*STDOUT, -format => 'Fasta');
}


if($opt_extermityOnly){
  print "option extremities activated => good for "
}
if( ($opt_upstreamRegion or $opt_downRegion)  and ! $opt_extermityOnly){print "I don't use upstream region option without the \"extermity only\" option (ext) to avoid issue with multi-features !"; exit;}

print "We will extract the $opt_type sequences.\n";
$opt_type=lc($opt_type);

if($codonTable<0 and $codonTable>25){
  print "$codonTable codon table is not a correct value. It should be between 0 and 25 (0,23 and 25 can be problematic !)\n";
}


##### MAIN ####
#### read gff file and save info in memory
######################
### Parse GFF input #
print "Reading file $opt_gfffile\n";
my ($hash_omniscient, $hash_mRNAGeneLink) = BILS::Handler::GFF3handler->slurp_gff3_file_JD($opt_gfffile);
print "Parsing Finished\n";
### END Parse GFF input #
#########################

my $hash_l1_grouped = group_l1features_from_omniscient($hash_omniscient);

#### read fasta
my $nbFastaSeq=0;
my $db = Bio::DB::Fasta->new($opt_fastafile);
print ("Genome fasta parsed\n");

foreach my $seqname (keys %{$hash_l1_grouped}) {
  foreach my $feature_l1 (@{$hash_l1_grouped->{$seqname}}) {
    my $id_l1=lc($feature_l1->_tag_value('ID'));
    my $name=undef;

    if ($feature_l1->has_tag('Name')){
      $name = $feature_l1->_tag_value('Name');
    }
    elsif($feature_l1->has_tag('gene')){
      $name = $feature_l1->_tag_value('gene');
    }

    if($opt_type eq lc($feature_l1->primary_tag()) ){

      #Handle Header
      my $header=$id_l1;
      if($name){
        $header.="|Name=".$name;
      }
      $seqname =~ tr/|/_/;
      $header.="|Seq_id=".$seqname."|type=".$opt_type;

      my @ListSeq=($feature_l1);
      my ($seqObj, $info) = extract_sequence(\@ListSeq, $db, $opt_extermityOnly, $opt_upstreamRegion, $opt_downRegion);
      if($info){
        $header.="|".$info;
      }
      $seqObj->id($header);

      print_seqObj($ostream, $seqObj, $opt_AA, $codonTable);
    }

    #################
    # == LEVEL 2 == #
    #################
    foreach my $ptag_l2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
         
      if ( exists ($hash_omniscient->{'level2'}{$ptag_l2}{$id_l1} ) ){
        foreach my $feature_l2 ( @{$hash_omniscient->{'level2'}{$ptag_l2}{$id_l1}}) {

          #For Header
          my $id_l2  = lc($feature_l2->_tag_value('ID'));
          if ($feature_l2->has_tag('Name') and ! $name){
            $name = $feature_l2->_tag_value('Name');
          }
          elsif($feature_l2->has_tag('gene') and ! $name){
            $name = $feature_l2->_tag_value('gene');
          }

          #Handle Header
          my $header=$id_l2."|".$id_l1;
          if($name){
            $header.="|Name=".$name;
          }
          $seqname =~ tr/|/_/;
          $header.="|Seq_id=".$seqname."|type=".$opt_type;

          if($opt_type eq $ptag_l2){
            my @ListSeq=($feature_l2);
            my ($seqObj, $info) = extract_sequence(\@ListSeq, $db, $opt_extermityOnly, $opt_upstreamRegion, $opt_downRegion);
            if($info){
              $header.="|".$info;
            }
            $seqObj->id($header);
            
            print_seqObj($ostream, $seqObj, $opt_AA, $codonTable);
          }

          #################
          # == LEVEL 3 == #
          #################
          foreach my $ptag_l3 (keys %{$hash_omniscient->{'level3'}}){
            if ( exists ($hash_omniscient->{'level3'}{$ptag_l3}{$id_l2} ) ){
              
              if( $opt_type eq $ptag_l3 ){
                my ($seqObj, $info) = extract_sequence(\@{$hash_omniscient->{'level3'}{$ptag_l3}{$id_l2}}, $db, $opt_extermityOnly, $opt_upstreamRegion, $opt_downRegion);
                if($info){
                  $header.="|".$info;
                }
                $seqObj->id($header);
                #print 
                print_seqObj($ostream, $seqObj, $opt_AA, $codonTable);
              }
            }
          }
        }
      }
    }
  }
}

#END
print "usage: $0 @copyARGV\n";

if($opt_upstreamRegion and $opt_downRegion){
  print "$nbFastaSeq $opt_type converted in fasta with $opt_upstreamRegion upstream nucleotides and $opt_downRegion downstream nucleotides.\n";
}
elsif($opt_upstreamRegion){
  print "$nbFastaSeq $opt_type converted in fasta with $opt_upstreamRegion upstream nucleotides.\n";
}
elsif($opt_downRegion){
  print "$nbFastaSeq $opt_type converted in fasta with $opt_downRegion downstream nucleotides.\n";
}
else{
  print "$nbFastaSeq $opt_type converted in fasta.\n";
}

my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job done in $run_time seconds\n";

#######################################################################################################################
        ####################
         #     METHODS    #
          ################
           ##############
            ############
             ##########
              ########
               ######
                ####
                 ##          

sub extract_sequence{
  my($feature_list, $db, $extermityOnly, $opt_upstreamRegion, $opt_downRegion)=@_;

  my @sortedList = sort {$a->start <=> $b->start} @$feature_list;
  my $sequence="";
  my $info="";
  if($extermityOnly){

    my $start = $sortedList[0]->start;
    my $end = $sortedList[$#sortedList]->end;

    #5'
    if($opt_upstreamRegion){
      #negative strand
      if($sortedList[0]->strand eq "-1" or $sortedList[0]->strand eq "-"){
        $end=$end+$opt_upstreamRegion;
        
        #get info
        if($end > $db->length($sortedList[0]->seq_id) ){
          $info.="5'extra=".($db->length($sortedList[0]->seq_id)-$end-$opt_upstreamRegion)."nt" ;
        }
        else{$info.="5'extra=".$opt_upstreamRegion."nt";}
      }
      else{
        $start=$start-$opt_upstreamRegion;        

        if($start < 0){
          $info.="5'extra=".($start+$opt_upstreamRegion)."nt";
        }
        else{
          $info.="5'extra=".$opt_upstreamRegion."nt";
        }
      }
    }
    #3'
    if($opt_downRegion){
      if($sortedList[0]->strand eq "-1" or $sortedList[0]->strand eq "-"){   
        $start=$start-$opt_upstreamRegion;

        #get info
        if($start < 0){
          $info.="3'extra=".($start+$opt_upstreamRegion)."nt";
        }
        else{$info.="3'extra=".$opt_upstreamRegion."nt";}
      }
      else{
        $end=$end+$opt_upstreamRegion;

        #get info
        if($end > $db->length($sortedList[0]->seq_id) ){
          $info.="3'extra=".$db->length($sortedList[0]->seq_id)-$end-$opt_upstreamRegion."nt" ;
        }
        else{$info.="3'extra=".$opt_upstreamRegion."nt";}
      }
    }

    $sequence = $db->subseq($sortedList[0]->seq_id, $start, $end);
  }
  else{
      foreach my $feature ( @sortedList ){
      $sequence .= $db->subseq($feature->seq_id,$feature->start,$feature->end);
    }
  }

  #create sequence object
  my $seq  = Bio::Seq->new( '-format' => 'fasta' , -seq => $sequence);
  
  #check if need to be reverse complement
  if($sortedList[0]->strand eq "-1" or $sortedList[0]->strand eq "-"){
    $seq=$seq->revcom;
  }

  return $seq,$info ;
}

sub print_seqObj{
  my($ostream, $seqObj, $opt_AA, $codonTable) = @_;

  $nbFastaSeq++;
  
  if($opt_AA){ #translate if asked
      my $transObj = $seqObj->translate(-CODONTABLE_ID => $codonTable);
      $ostream->write_seq($transObj);  
   }
  else{
    $ostream->write_seq($seqObj);                
  }
}


__END__

=head1 NAME

gff3_extract_cds_sequences.pl -
This script extract sequence in feasta format from gff file. You can extract the fasta of any kind of feature define by the 2th column in the gff file.
The result is written to the specified output file, or to STDOUT.

The Header are formated like that:
>mRNA_ID|gene_ID|Name=NAME|Seq_id=Chromosome_ID|type=cds|5'extra=VALUE

/!\ mRNA_ID not displayed when extracting gene. 
Name is optional and will be written only if the Name attribute exists in th gff.
type will be the feature type extracted. 
5'extra or 3'extra is otpional, according to the use of the upstream and downstream option.

=head1 SYNOPSIS

    ./gff3_extract_cds_sequences.pl -g=infile.gff -f=infile.fasta  [ -o outfile ]
    ./gff3_extract_cds_sequences.pl --help

=head1 OPTIONS

=over 8

=item B<-g>, B<--gff> or B<-ref>

Input GFF3 file that will be read (and sorted)

=item B<-f> or B<--fasta> 

Input fasta file.

=item B<-t> 

Define the feature you want to extract the sequnce from. By deafault it's 'cds'. Most common choice are: gene,mrna,exon,cds,trna,three_prime_utr,five_prime_utr.

=item B<-p>, B<--protein> or B<--aa>

Will translate the extracted sequence in Amino acid. By default the codon table used is the 1 (Standard). See codon table option for more options.

=item B<--codon>, B<--table> or B<--ct>

Allow to choose another type of codon table for the translation.

=item B<-e> or B<--ext>

This option called "extremities" allows dealing with multifeature like cds or exon, to extract the full sequence from start extremity to the end extremity, i.e with introns.
Use of that option with exon will give the same result as extract the mrna sequence.
Use of that option on cds will give the cdna wihtout the untraslated sequences.

=item B<-u>, B<--up>, B<-5>, B<--five> or B<-upstream>

Integer. It will take that number of nucleotide in more at the 5' extremity. Option "e" must be activated to use this option (Why ? to avoid to extract intronic/overlaping sequence in case of feature spread over several locations (exon,cds,utrs)).

=item B<-d>, B<--do>, B<-3>, B<--three>, B<-down> or B<-downstream>

Integer. It will take that number of nucleotide in more at the 3' extremity. Option "e" must be activated to use this option (Why ? to avoid to extract intronic/overlaping sequence in case of feature spread over several locations (exon,cds,utrs)).

=item B<--cdna>

This extract the cdna sequence (i.e transcribed sequence (devoid of introns, but containing untranslated exons)). It correspond to extract the exons sequences.

=item B<-o> or B<--output>

Output GFF file.  If no output file is specified, the output will be
written to STDOUT.

=item B<-h> or B<--help>

Display this helpful text.

=back

=cut
