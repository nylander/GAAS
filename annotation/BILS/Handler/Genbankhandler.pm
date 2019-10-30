#!/usr/bin/perl -w

package BILS::Handler::Genbankhandler ;

use strict;
use Clone 'clone';

use Exporter qw(import);

our $VERSION     = 1.00;
our @ISA         = qw(Exporter);
our @EXPORT_OK   = qw();
our %EXPORT_TAGS = ( DEFAULT => [qw()],
                 Ok    => [qw()]);
=head1 SYNOPSIS



=head1 DESCRIPTION

	A library about genbank format
	Inherits from

=cut

##################
#
#	About Genbank Contants
#
##################


use constant GENBANK_DIVISION => {
    PRI => 'primate sequences', ROD => 'rodent sequences', MAM => 'other mammalian sequences', VRT => 'other vertebrate sequences', INV => 'invertebrate sequences', PLN => 'plant, fungal, and algal sequences',
    BCT => 'bacterial sequences', VRL => 'viral sequences', PHG => 'bacteriophage sequences', SYN => 'synthetic sequences', UNA => 'unannotated sequences', EST => 'EST sequences (expressed sequence tags)',
    PAT => 'patent sequences', STS => 'STS sequences (sequence tagged sites)', GSS => 'GSS sequences (genome survey sequences)', HTG => 'HTG sequences (high-throughput genomic sequences)',
    HTC => 'unfinished high-throughput cDNA sequencing', ENV => 'environmental sampling sequences',
    };



# Save in omniscient hash (sorted in a specific way (3 levels)) a whole gff3 file
sub slurp_genbank_file_JD {

	my ($self, $file) = @_  ;

	my $gtfio = Bio::Tools::GFF->new(-file => $file, -gff_version => 3);

	my %mRNAGeneLink;
	my %omniscient;

	while( my $feature = $gtfio->next_feature()) {
		manage_one_feature($feature, \%omniscient, \%mRNAGeneLink);
    }

	return \%omniscient, \%mRNAGeneLink	;
}

sub create_omniscient_from_feature_list {
	my ($ref_featureList)=@_;

	my %mRNAGeneLink;
	my %omniscient;

	foreach my $feature (@{$ref_featureList}) {
		manage_one_feature($feature, \%omniscient, \%mRNAGeneLink);
    }

	return \%omniscient, \%mRNAGeneLink	;
}

sub manage_one_feature{

	my ($feature, $omniscient, $mRNAGeneLink)=@_;

		my $seq_id = $feature->seq_id;					#col1
		my $source_tag = $feature->source_tag;			#col2
		my $primary_tag = lc($feature->primary_tag);	#col3
		my $start = $feature->start;					#col4
		my $end = $feature->end;						#col5
		my $score = $feature->score;					#col6
		my $strand = $feature->strand;					#col7
		my $frame = $feature->frame;					#col8
		#Attribute-value								#col9
		my $id= undef;
		my $parent= undef;

		####################################################
	    ########## Manage feature WITHOUT parent ###########	== LEVEL1 ==
	    ####################################################
	    if ($primary_tag eq "gene" ) {
	    	#get ID
	    	my @values = $feature->get_tag_values('ID');
	    	if($#values == -1){print "error !! No ID attribute found for the feature $feature->gff_string()\n";}
			$id = lc(shift @values) ;

	        $omniscient->{"level1"}{$primary_tag}{$id}=$feature;
	        next();
	    }

      	###################################################
      	########## Manage feature WITHout CHILD  ##########		== LEVEL3 ==
      	###################################################

      	elsif ( ($primary_tag eq "cds") or ($primary_tag eq "exon") or ($primary_tag eq "stop_codon") or ($primary_tag eq "start_codon") or ($primary_tag eq "three_prime_utr") or ($primary_tag eq "five_prime_utr") or ($primary_tag eq "utr"),
      	 or ($primary_tag eq "selenocysteine") or ($primary_tag eq "non_canonical_three_prime_splice_site") or ($primary_tag eq "non_canonical_five_prime_splice_site") or ($primary_tag eq "stop_codon_read_through") ){
      		my @values = $feature->get_tag_values('Parent');
	    	if($#values == -1){print "error !! No Parent attribute found for the feature $feature->gff_string()\n";}
			$parent = lc(shift @values) ;

      		push (@{$omniscient->{"level3"}{$primary_tag}{$parent}}, $feature);
      	}
      	##############################################
      	########## Manage feature the rest  ##########		== LEVEL2 ==
      	##############################################
      	elsif( ($primary_tag eq "mrna") or ($primary_tag eq "ncrna") or ($primary_tag eq "mirna") or ($primary_tag eq "lcrna") or ($primary_tag eq "rrna") or ($primary_tag eq "srp_rna"),
      		or ($primary_tag eq "snrna") or ($primary_tag eq "lincrna") or ($primary_tag eq "trna") or ($primary_tag eq "trna_pseudogene") or ($primary_tag eq "snorna") or ($primary_tag eq "misc_rna"),
      		or ($primary_tag eq "rnase_p_rna") ){
    		#get ID
	    	my @values = $feature->get_tag_values('ID');
	    	if($#values == -1){print "error2 !! No ID found for the feature $feature->gff_string()\n";}
			$id = lc(shift @values) ;

			my @values = $feature->get_tag_values('Parent');
	    	if($#values == -1){print "error2 !! No Parent attribute found for the feature $feature->gff_string()\n";}
			$parent = lc(shift @values);

  			if (! exists ($mRNAGeneLink->{$id})){ # keep track of link between level2->leve1
				$mRNAGeneLink->{$id}=$parent;
	 		}

      		push (@{$omniscient->{"level2"}{$primary_tag}{$parent}}, $feature);
      	}

      	else{
      		print "gff reader : $primary_tag still not taken in account ! Please modify the code to define on of the three level of this feature.\n";
      	}

}



1;
