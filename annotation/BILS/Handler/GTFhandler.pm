#!/usr/bin/perl -w

package BILS::Handler::GTFhandler ;

use Data::Dumper;
use strict;
use BILS::Handler::GXFhandler qw(:Ok);
use Exporter qw(import);

our $VERSION     = 1.00;
our @ISA         = qw(Exporter);
our @EXPORT_OK   = qw( manage_level3_attributes_ID_parent sort_gtf_feature_list);
our %EXPORT_TAGS = ( DEFAULT => [qw()],
                 Ok    => [qw( manage_level3_attributes_ID_parent sort_gtf_feature_list)]);
=head1 SYNOPSIS



=head1 DESCRIPTION

	A library to convert 
	Inherits from 

	Dont take in account repeat and multi parent feature!!!
	
=cut	

# Save in omniscient hash (sorted in a specific way (3 levels)) a whole gtf file
sub slurp_gtf_file_JD {
	
	my ($self, $file) = @_  ;
	
	my $gtfio = Bio::Tools::GFF->new(-file => $file, -gff_version => 2.5);	

	my %mRNAGeneLink;
	my %omniscient;

	while( my $feature = $gtfio->next_feature()) {
		manage_one_feature($feature, \%omniscient, \%mRNAGeneLink);
    }
    
	return \%omniscient, \%mRNAGeneLink	;
}

sub sort_gtf_feature_list {
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
		my $primary_tag = $feature->primary_tag;		#col3
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
	    if (lc($primary_tag) eq "gene" ) {
	    	#get ID
	    	my @values = $feature->get_tag_values('gene_id');
	    	if($#values == -1){print "error !! No gene_id found for the feature $feature->gff_string()\n";}
			$id = shift @values ;
    
	        $omniscient->{"level1"}{lc($primary_tag)}{lc($id)}=$feature;
	        next();
	    }

      	###################################################
      	########## Manage feature WITHout CHILD  ##########		== LEVEL3 ==
      	###################################################

      	elsif ( (lc($primary_tag) eq "cds" ) or (lc($primary_tag) eq "exon" ) or (lc($primary_tag) eq "stop_codon" ) or (lc($primary_tag) eq "start_codon" ) or (lc($primary_tag) eq "utr" ) or (lc($primary_tag) eq "selenocysteine" ) ){
      		if($feature->has_tag('transcript_id')){
	      		$parent = $feature->_tag_value('transcript_id');
	      		push (@{$omniscient->{"level3"}{lc($primary_tag)}{lc($parent)}}, $feature);
	      	}
	      	else{# Assume there is only one isoform !! We take transcript_id from annother level3 feature already registered of the same feature group
	      		warn "GTFhandler Warning !! transcript_id was not found for the following feature: \n".$feature->gff_string."\n";
	      		print "We assume only one isoform exists and will try to get the transcript_id from annother level3 feature having the same level1 gene_id\n";
	      		
	      		if($feature->has_tag('gene_id')){ # Check if gene id exists
		      		my $gene_id = $feature->_tag_value('gene_id');
		      		print "My geneid=".$gene_id."\n";
			      	foreach my $idm (keys %{$omniscient->{"level3"}{'exon'}}){
			      		foreach my $f ( @{$omniscient->{"level3"}{'exon'}{$idm}}){
			      			
			      			if($f->has_tag('gene_id')){ # Check if gene id exists
				      			my $f_id=$f->_tag_value('gene_id');
				      			
				      			if($f_id eq $gene_id){

				      				if($f->has_tag('transcript_id')){
				      					$parent = $f->_tag_value('transcript_id');
				      					create_or_replace_tag($feature,'transcript_id',$parent);
				      					push (@{$omniscient->{"level3"}{lc($primary_tag)}{lc($parent)}}, $feature);
				      					print "now: ".$feature->gff_string."\n\n";
				      					last;
				      				}
								}
							}				
		      			}
		      			if($parent){last;}
		      		}
		      	}
		      	else{ #No geneid we cannot look for the thranscript ID
		      		print "No luck ! No gene id neither for this feature. A last try will be to get a protein_id, hoping that it will be similar as the transcript_id (i.e JGI may do that)\n";
		      		if($feature->has_tag('protein_id')){ # Check if protein_id exists
		      			$parent = $feature->_tag_value('protein_id');
		      			create_or_replace_tag($feature,'transcript_id',$parent);
		      			print "now: ".$feature->gff_string."\n\n";
		      			push (@{$omniscient->{"level3"}{lc($primary_tag)}{lc($parent)}}, $feature);
		      		}
		      		else{
		      			print "Don't push me too much, you bother me !! I'm really nice so I will try a last trick. I will check if a name attribute exists and if it is similar to another feature from the same level I will check if this other feature has what I'm looking for.\n";
		      			if($feature->has_tag('name')){ # Check if gene id exists
				      		my $gene_name = $feature->_tag_value('name');
				      		print "name used =".$gene_name."\n";
					      	foreach my $idm (keys %{$omniscient->{"level3"}{'exon'}}){
					      		foreach my $f ( @{$omniscient->{"level3"}{'exon'}{$idm}}){
					      			
					      			if($f->has_tag('name')){ # Check if gene id exists
						      			my $f_name=$f->_tag_value('name');
						      			
						      			if($f_name eq $gene_name){

						      				if($f->has_tag('transcript_id')){
						      					$parent = $f->_tag_value('transcript_id');
						      					create_or_replace_tag($feature,'transcript_id',$parent);
						      					print "now: ".$feature->gff_string."\n\n";
						      					push (@{$omniscient->{"level3"}{lc($primary_tag)}{lc($parent)}}, $feature);
						      					last;
						      				}
						      				elsif($feature->has_tag('protein_id')){ # Check if protein_id exists
		      									$parent = $feature->_tag_value('protein_id');
		      									create_or_replace_tag($feature,'transcript_id',$parent);
		      									push (@{$omniscient->{"level3"}{lc($primary_tag)}{lc($parent)}}, $feature);
		      									print "now: ".$feature->gff_string."\n\n";
		      									last;
		      								}	
										}
									}				
				      			}
				      			if($parent){last;}
				      		}
		      			}
		      			else{
			      			print "Sorry, we cannot manage this feature. We tried everything we can.\n"; exit;
			      		}
		      		}
		      	}
	      	}
      	}
      	##############################################
      	########## Manage feature the rest  ##########		== LEVEL2 ==
      	##############################################
      	elsif( (lc($primary_tag) eq "mrna" ) or (lc($primary_tag) eq "ncrna" ) or (lc($primary_tag) eq "mirna" ) or (lc($primary_tag) eq "lcrna" ) or (lc($primary_tag) eq "transcript" ) ){
    		#get ID
    		if($feature->has_tag('transcript_id')){ # Check if transcript id exists
	    		$id = $feature->_tag_value('transcript_id');
	    	}else{warn "GTFhandler Warning !! transcript_id was not found for the following feature: \n".$feature->gff_string."\n";}

	    	if($feature->has_tag('gene_id')){ # Check if gene id exists
				$parent = $feature->_tag_value('gene_id');
			}else{warn "GTFhandler Warning !! transcript_id was not found for the following feature: \n".$feature->gff_string."\n";}

  			if (! exists ($mRNAGeneLink->{lc($id)})){ # keep track of link between level2->leve1
				$mRNAGeneLink->{lc($id)}=lc($parent);
	 		}
      		push (@{$omniscient->{"level2"}{lc($primary_tag)}{lc($parent)}}, $feature);
      	}
      	else{
      		print "gtf reader : $primary_tag still not taken in account ! Please modify the code to define this feature as one of the three level possible.\n";
      	}	

}

##==============================================================

sub manage_level3_attributes_ID_parent {
	my $features = \@_;

	# Take all features and bin them by transcript_id	
	my %bin = _group_features_by_transcript_and_seq_id(@$features);

	my @answer = ();

	# Some annotations (like lift-overs) can be split across scaffolds/contigs
	my $is_split = 0;

	# Check wether this gene is split across scaffolds
	if ( (keys %bin) > 1) {
		$is_split = 1;
	}		

	my $seq_counter = 0;

	while (my ($seq_id,$transcript_hash) = each %bin) {

		# All other features to array, neatly ordered and modified where needed
		while ( my ($transcript_id, $values) = each(%$transcript_hash) ) {

			# Features other than mRNAs and genes need to have new stable IDs, let's build some...
			my $feature_counter = 0;
			my $primary_tag_before;
			foreach my $f (@$values) {
				
				#manage counter of feature type
				my $primary_tag = lc($f->primary_tag);
				if ($primary_tag_before ne $primary_tag){$feature_counter = 0;} # reinitialise the counter if as example : we had exon and now we have cds (value previously sorted)
 				$feature_counter += 1;

				# Manage case where features spread over several seq
				my $appendix ="";
				if ($is_split == 1) {
					my $seqid=$f->seq_id();
					$appendix = "_partial_part-$seqid" ;   ############# //// !!!! \\\\\ _partial_part- is used in kraken_statMap script to the total lenght of the original mRMA before to be split.
				}
				my $correct_transcript_id="$transcript_id$appendix";

				#manage id name
 				my $id=undef;			
 				if ($primary_tag eq "exon"){
 					if ($f->has_tag('exon_id')){
 				 		 my $id  = $f->_tag_value('exon_id');
 				 	}
 				}
 				elsif($primary_tag eq "cds"){
 				 	if ($f->has_tag('protein_id')){
 				 		my $id = $f->_tag_value('protein_id'). "-" . $f->primary_tag . "-" . $feature_counter ;
 				 	}
 				 }
 				if($id eq undef){
 					$id = $correct_transcript_id . "-" . $f->primary_tag . "-" . $feature_counter ; # A crummy new stable ID...
 				}


 				create_or_replace_tag($f, 'ID', $id );
 				create_or_replace_tag($f, 'Parent', $correct_transcript_id );

 				push(@answer,$f);

 				$primary_tag_before=$primary_tag;
 			}
		}
 	}
	return \@answer;
}

1;
