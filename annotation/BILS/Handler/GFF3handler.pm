#!/usr/bin/perl -w

package BILS::Handler::GFF3handler ;


use strict;
use warnings;
use Data::Dumper;
use Exporter qw(import);
use Sort::Naturally;
use URI::Escape;
use Bio::Seq;

our $VERSION     = 1.00;
our @ISA         = qw(Exporter);
our @EXPORT_OK   = qw(check_record_positions complement_omniscients_by_size convert_omniscient_to_ensembl_style remove_l2_related_feature l2_identical group_l1IDs_from_omniscient complement_omniscients rename_ID_existing_in_omniscient keep_only_uniq_from_list2 check_gene_overlap_at_CDSthenEXON location_overlap_update location_overlap print_omniscient_as_match nb_feature_level1 check_gene_positions find_overlap_between_geneFeature_and_sortBySeqId sort_by_seq_id sort_by_seq webapollo_compliant extract_cds_sequence group_l1features_from_omniscient create_omniscient_from_idlevel2list get_feature_l2_from_id_l2_l1 remove_omniscient_elements_from_level2_feature_list remove_omniscient_elements_from_level2_ID_list featuresList_identik group_features_from_omniscient featuresList_overlap check_level1_positions check_level2_positions info_omniscient fil_cds_frame exists_keys remove_element_from_omniscient append_omniscient merge_omniscients remove_omniscient_elements_from_level1_id_list fill_omniscient_from_other_omniscient_level1_id print_omniscient_from_level1_id_list subsample_omniscient_from_level1_id_list check_if_feature_overlap remove_tuple_from_omniscient print_ref_list_feature print_omniscient create_or_replace_tag remove_element_from_omniscient_attributeValueBased get_longest_cds_level2);
our %EXPORT_TAGS = ( DEFAULT => [qw()],
                 Ok    => [qw(check_record_positions complement_omniscients_by_size convert_omniscient_to_ensembl_style remove_l2_related_feature l2_identical group_l1IDs_from_omniscient complement_omniscients rename_ID_existing_in_omniscient keep_only_uniq_from_list2 check_gene_overlap_at_CDSthenEXON location_overlap_update location_overlap print_omniscient_as_match nb_feature_level1 check_gene_positions find_overlap_between_geneFeature_and_sortBySeqId sort_by_seq_id sort_by_seq webapollo_compliant extract_cds_sequence group_l1features_from_omniscient create_omniscient_from_idlevel2list get_feature_l2_from_id_l2_l1 remove_omniscient_elements_from_level2_feature_list remove_omniscient_elements_from_level2_ID_list featuresList_identik group_features_from_omniscient featuresList_overlap check_level1_positions check_level2_positions info_omniscient fil_cds_frame exists_keys remove_element_from_omniscient append_omniscient merge_omniscients remove_omniscient_elements_from_level1_id_list fill_omniscient_from_other_omniscient_level1_id print_omniscient_from_level1_id_list subsample_omniscient_from_level1_id_list check_if_feature_overlap remove_tuple_from_omniscient print_ref_list_feature print_omniscient create_or_replace_tag remove_element_from_omniscient_attributeValueBased get_longest_cds_level2)]);
=head1 SYNOPSIS



=head1 DESCRIPTION

	A library to convert
	Inherits from

	Dont take in account repeat and multi parent feature!!!

=cut

#				   +------------------------------------------------------+
#				   |+----------------------------------------------------+|
#				   || 					Print Methods 					 ||
#				   |+----------------------------------------------------+|
#				   +------------------------------------------------------+

# 1) Original print
###################
sub print_ref_list_feature {

	my ($list, $gffout) = @_  ;

	foreach my $feature (@$list) {
		$gffout->write_feature($feature);
	}
}

# omniscient is a hash containing a whole gXf file in memory sorted in a specific way (3 levels)
sub print_omniscient{

	my ($hash_omniscient, $gffout) = @_  ;

	#uri_decode_omniscient($hash_omniscient);

### OLD FASHION GOING TRHOUGH LEVEL1
	#foreach my $primary_tag_l1 ( sort {$a <=> $b or $a cmp $b} keys %{$hash_omniscient->{'level1'}}){ # primary_tag_l1 = gene or repeat etc...
	#	foreach my $id_tag_key_level1 ( sort { $hash_omniscient->{'level1'}{$primary_tag_l1}{$a}->start <=> $hash_omniscient->{'level1'}{$primary_tag_l1}{$b}->start } keys %{$hash_omniscient->{'level1'}{$primary_tag_l1}} ) { #sort by position

### NEW FASHION GOING TRHOUGH LEVEL1 - Have to first create a hash of seq_id -> level1_feature , then we can go through in alphanumerical order.
	# sort by seq id
	my %hash_sortBySeq;
	foreach my $tag_level1 ( keys %{$hash_omniscient->{'level1'}}){
	  foreach my $level1_id ( keys %{$hash_omniscient->{'level1'}{$tag_level1}}){
	    my $position=$hash_omniscient->{'level1'}{$tag_level1}{$level1_id}->seq_id;
	    push (@{$hash_sortBySeq{$position}{$tag_level1}}, $hash_omniscient->{'level1'}{$tag_level1}{$level1_id});
	  }
	}

	# Read by seqId to sort properly the output by seq ID
	# sort { (($a =~ /(\d+)$/)[0] || 0) <=> (($b =~ /(\d+)$/)[0] || 0) will provide sorting liek that: contig contig1 contig2 contig3 contig10 contig11 contig22 contig100 contig101
	foreach my $seqid (sort { (($a =~ /(\d+)$/)[0] || 0) <=> (($b =~ /(\d+)$/)[0] || 0) } keys %hash_sortBySeq){ # loop over all the feature level1

	  	foreach my $primary_tag_l1 (sort {$a cmp $b} keys %{$hash_sortBySeq{$seqid}}){

		    foreach my $feature_l1 ( sort { ncmp ($a->start.$a->end.$a->_tag_value('ID'), $b->start.$b->end.$b->_tag_value('ID') ) } @{$hash_sortBySeq{$seqid}{$primary_tag_l1}}){
			    my $id_tag_key_level1 = lc($feature_l1->_tag_value('ID'));
				$gffout->write_feature($hash_omniscient->{'level1'}{$primary_tag_l1}{$id_tag_key_level1}); # print feature

				#################
				# == LEVEL 2 == #
				#################
				foreach my $primary_tag_l2 (sort {$a cmp $b} keys %{$hash_omniscient->{'level2'}}){ # primary_tag_l2 = mrna or mirna or ncrna or trna etc...

					if ( exists_keys( $hash_omniscient, ('level2', $primary_tag_l2, $id_tag_key_level1) ) ){
						foreach my $feature_level2 ( sort { ncmp ($a->start.$a->end.$a->_tag_value('ID'), $b->start.$b->end.$b->_tag_value('ID') ) } @{$hash_omniscient->{'level2'}{$primary_tag_l2}{$id_tag_key_level1}}) {
							$gffout->write_feature($feature_level2);

							#################
							# == LEVEL 3 == #
							#################
							my $level2_ID = lc($feature_level2->_tag_value('ID'));

							######
							# FIRST EXON
							if ( exists_keys( $hash_omniscient, ('level3', 'exon', $level2_ID) ) ){
								foreach my $feature_level3 ( sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{'exon'}{$level2_ID}}) {
									$gffout->write_feature($feature_level3);
								}
							}
							###########
							# SECOND CDS
							if ( exists_keys( $hash_omniscient, ('level3', 'cds', $level2_ID) ) ){
								foreach my $feature_level3 ( sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{'cds'}{$level2_ID}}) {
									$gffout->write_feature($feature_level3);
								}
							}

							############
							# THEN ALL THE REST
							foreach my $primary_tag_l3 (sort {$a cmp $b} keys %{$hash_omniscient->{'level3'}}){ # primary_tag_l3 = cds or exon or start_codon or utr etc...
								if (($primary_tag_l3 ne 'cds') and ($primary_tag_l3 ne 'exon')) {
									if ( exists_keys( $hash_omniscient, ('level3', $primary_tag_l3, $level2_ID) ) ){
										foreach my $feature_level3 ( sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{$primary_tag_l3}{$level2_ID}}) {
											$gffout->write_feature($feature_level3);
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

# omniscient is a hash containing a whole gXf file in memory sorted in a specific way (3 levels)
sub print_omniscient_as_match{

	my ($hash_omniscient, $gffout) = @_  ;

	#uri_decode_omniscient($hash_omniscient);

### OLD FASHION GOING TRHOUGH LEVEL1
#	foreach my $primary_tag_l1 ( sort {$a <=> $b or $a cmp $b} keys %{$hash_omniscient->{'level1'}}){ # primary_tag_l1 = gene or repeat etc...
#		foreach my $id_tag_key_level1 ( sort { $hash_omniscient->{'level1'}{$primary_tag_l1}{$a}->start <=> $hash_omniscient->{'level1'}{$primary_tag_l1}{$b}->start } keys %{$hash_omniscient->{'level1'}{$primary_tag_l1}} ) { #sort by position

### NEW FASHION GOING TRHOUGH LEVEL1 - Have to first create a hash of seq_id -> level1_feature , then we can go tthrough in alphanumerical order.
	# sort by seq id
	my %hash_sortBySeq;
	foreach my $tag_level1 ( keys %{$hash_omniscient->{'level1'}}){
	  foreach my $level1_id ( keys %{$hash_omniscient->{'level1'}{$tag_level1}}){
	    my $position=$hash_omniscient->{'level1'}{$tag_level1}{$level1_id}->seq_id;
	    push (@{$hash_sortBySeq{$position}{$tag_level1}}, $hash_omniscient->{'level1'}{$tag_level1}{$level1_id});
	  }
	}

	#Read by seqId to sort properly the output by seq ID
	foreach my $seqid ( sort { (($a =~ /(\d+)$/)[0] || 0) <=> (($b =~ /(\d+)$/)[0] || 0) } keys %hash_sortBySeq){ # loop over all the feature level1
	  	foreach my $primary_tag_l1 (sort {$a cmp $b} keys %{$hash_sortBySeq{$seqid}}){

	#################
	# == LEVEL 1 == #
	#################
		    foreach my $feature_l1 ( sort {$a->start <=> $b->start} @{$hash_sortBySeq{$seqid}{$primary_tag_l1}}){
			    my $id_tag_key_level1 = lc($feature_l1->_tag_value('ID'));

				if($primary_tag_l1 =~ "match"){
					$gffout->write_feature($hash_omniscient->{'level1'}{$primary_tag_l1}{$id_tag_key_level1}); # print feature
				}
				#################
				# == LEVEL 2 == #
				#################
				foreach my $primary_tag_l2 (sort {$a cmp $b} keys %{$hash_omniscient->{'level2'}}){ # primary_tag_l2 = mrna or mirna or ncrna or trna etc...

					if ( exists_keys( $hash_omniscient, ('level2', $primary_tag_l2, $id_tag_key_level1) ) ){
						foreach my $feature_level2 ( sort {$a->start <=> $b->start} @{$hash_omniscient->{'level2'}{$primary_tag_l2}{$id_tag_key_level1}}) {
							
							if($primary_tag_l2 =~ "match"){
								$gffout->write_feature($feature_level2);
							}	
							else{
								$feature_level2->primary_tag('match');
								if( $feature_level2->has_tag('Parent')){
									$feature_level2->remove_tag('Parent');
								}

								$gffout->write_feature($feature_level2);

								#################
								# == LEVEL 3 == #
								#################
								my $level2_ID = $feature_level2->_tag_value('ID');

								######
								# EXON
								if ( exists_keys( $hash_omniscient, ('level3', 'exon', lc($level2_ID)) ) ){
									my $current_start=0;
									foreach my $feature_level3 ( sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{'exon'}{lc($level2_ID)}}) {
										
										$current_start++;
										my $end=($feature_level3->end - $feature_level3->start)+$current_start;
										$feature_level3->primary_tag('match_part');
										
										if(! $feature_level3->has_tag('Target')){
											my @target=();
											create_or_replace_tag($feature_level3, "Target", "$level2_ID $current_start $end +"); # Target has value has to be a list correctly formated
										}
										$current_start=$end;

										$gffout->write_feature($feature_level3);
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

# omniscient is a hash containing a whole gXf file in memory sorted in a specific way (3 levels)
sub print_omniscient_from_level1_id_list {

	my ($hash_omniscient, $level_id_list, $gffout) = @_  ;

	#################
	# == LEVEL 1 == #
	#################
	foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...

		foreach my $id_tag_key_level1_raw (@$level_id_list){
			my $id_tag_key_level1 = lc($id_tag_key_level1_raw);
			if(exists ($hash_omniscient->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1})){

				#uri_encode_one_feature($hash_omniscient->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1});

				$gffout->write_feature($hash_omniscient->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1}); # print feature

				#################
				# == LEVEL 2 == #
				#################
				foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

					if ( exists ($hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1} ) ){
						foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}) {

							#uri_encode_one_feature($feature_level2);

							$gffout->write_feature($feature_level2);

							#################
							# == LEVEL 3 == #
							#################
							my $level2_ID ;
							if($feature_level2->has_tag('ID')){
								$level2_ID = lc($feature_level2->_tag_value('ID'));
							}
							elsif($feature_level2->has_tag('transcript_id')){
								$level2_ID = lc( $feature_level2->_tag_value('transcript_id'));
							}
							else{
								warn "Cannot retrieve the parent feature of the following feature: ".gff_string($feature_level2);
							}

							###########
							# Before tss
							if ( exists_keys($hash_omniscient,('level3','tss',$level2_ID)) ){
								foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{'tss'}{$level2_ID}}) {

									#uri_encode_one_feature($feature_level3);

									$gffout->write_feature($feature_level3);
								}
							}

							######
							# FIRST EXON
							if ( exists_keys($hash_omniscient,('level3','exon',$level2_ID)) ){
								foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{'exon'}{$level2_ID}}) {

									#uri_encode_one_feature($feature_level3);

									$gffout->write_feature($feature_level3);
								}
							}
							###########
							# SECOND CDS
							if ( exists_keys($hash_omniscient,('level3','cds',$level2_ID)) ){
								foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{'cds'}{$level2_ID}}) {

									#uri_encode_one_feature($feature_level3);

									$gffout->write_feature($feature_level3);
								}
							}

							###########
							# Last tts
							if ( exists_keys($hash_omniscient,('level3','tts',$level2_ID)) ){
								foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{'tts'}{$level2_ID}}) {

									#uri_encode_one_feature($feature_level3);

									$gffout->write_feature($feature_level3);
								}
							}

							###########
							# The rest
							foreach my $primary_tag_key_level3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
								if( ($primary_tag_key_level3 ne 'cds') and ($primary_tag_key_level3 ne 'exon') and ($primary_tag_key_level3 ne 'tss') and ($primary_tag_key_level3 ne 'tts')){
									if ( exists ($hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID} ) ){
										foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID}}) {

											#uri_encode_one_feature($feature_level3);

											$gffout->write_feature($feature_level3);
										}
									}
								}
							}
						}
					}
				}
			}
		}

	}
}

###############################
# METHOD RELATED TO WEBAPOLLO #
###############################
use constant CWA_skip_feature => { "non_canonical_three_prime_splice_site" => 1 , "non_canonical_five_prime_splice_site" => 2};
#Transform omniscient data to be Webapollo compliant
sub webapollo_compliant {
		my ($hash_omniscient) = @_  ;

	#################
	# == LEVEL 1 == #
	#################
	foreach my $primary_tag_l1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_l1 = gene or repeat etc...
		if(exists (CWA_skip_feature->{$primary_tag_l1})){delete $hash_omniscient->{'level1'}{$primary_tag_l1}; next;}
		foreach my $id_l1 (keys %{$hash_omniscient->{'level1'}{$primary_tag_l1}}){
			webapollo_rendering_l1($hash_omniscient->{'level1'}{$primary_tag_l1}{$id_l1});
			#################
			# == LEVEL 2 == #
			#################
			foreach my $primary_tag_l2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_l2 = mrna or mirna or ncrna or trna etc...
				if(exists (CWA_skip_feature->{$primary_tag_l2})){delete $hash_omniscient->{'level2'}{$primary_tag_l2}; next;}
				if ( exists ($hash_omniscient->{'level2'}{$primary_tag_l2}{$id_l1} ) ){
					foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_l2}{$id_l1}}) {
						webapollo_rendering_l2($feature_level2);

						#################
						# == LEVEL 3 == #
						#################
						my $level2_ID  = lc($feature_level2->_tag_value('ID'));

						foreach my $primary_tag_l3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_l3 = cds or exon or start_codon or utr etc...
							if(exists (CWA_skip_feature->{$primary_tag_l3})){delete $hash_omniscient->{'level3'}{$primary_tag_l3}; next;}
							if ( exists ($hash_omniscient->{'level3'}{$primary_tag_l3}{$level2_ID} ) ){
								foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{$primary_tag_l3}{$level2_ID}}) {
									webapollo_rendering_l3($feature_level3);
								}
							}
						}
					}
				}
			}
		}
	}
}

#follow webapollo description for a correct visualisation of data 
sub webapollo_rendering_l1 {
	my ($feature)=@_;

	my @tags = $feature->get_all_tags();
	
	## check Name attribute	
	my @f = grep /^\Qname\E$/i, @tags; #look at tag that match the whole string to NAME (case insensitive)
	my $name_tag = $f[0];
	if(! $name_tag){
		my $value = $feature->_tag_value('ID');
		$feature->add_tag_value('Name', $value);
	}
}

	my ($feature)=@_;

#follow webapollo description for a correct visualisation of data 
sub webapollo_rendering_l2 {
	my ($feature)=@_;

	## check primary tag
	my $primary_tag = lc($feature->primary_tag);

	my %corrections = (
			mrna => 'mRNA',
		);
	if ( exists $corrections{$primary_tag}) {
		$feature->primary_tag( $corrections{$primary_tag});
	}

	my @tags = $feature->get_all_tags();

	## check product/description attribute
	#	if($feature->has_tag('product')){
	my @f = grep /\Qproduct\E/i, @tags;
	my $product_tag = $f[0];
	if($product_tag){
		my @values = $feature->get_tag_values($product_tag);
		$feature->add_tag_value('description', @values);
		$feature->remove_tag($product_tag);
	}
}

#follow webapollo description for a correct visualisation of data 
sub webapollo_rendering_l3 {
	my ($feature)=@_;

	## check primary tag
	my $primary_tag = lc($feature->primary_tag);

	my %corrections = (
			cds => 'CDS',
			exon => 'exon',
			three_prime_utr => 'three_prime_UTR',
			five_prime_utr => 'five_prime_UTR',
			utr => 'UTR',
		);
	if ( exists $corrections{$primary_tag}) {
		$feature->primary_tag( $corrections{$primary_tag});
	}

	my @tags = $feature->get_all_tags();

	foreach my $tag (@tags){
		if(lc($tag) ne 'id' and lc($tag) ne 'parent'){
			$feature->remove_tag($tag);
		}
	}
}

#Transform omniscient data to be embl compliant
sub embl_compliant {
		my ($hash_omniscient) = @_  ;

	#################
	# == LEVEL 1 == #
	#################
	foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
		foreach my $id_tag_key_level1 (keys %{$hash_omniscient->{'level1'}{$primary_tag_key_level1}}){
			embl_rendering($hash_omniscient->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1});

			#################
			# == LEVEL 2 == #
			#################
			foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

				if ( exists ($hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1} ) ){
					foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}) {
						embl_rendering($feature_level2);

						#################
						# == LEVEL 3 == #
						#################
						my $level2_ID  = lc($feature_level2->_tag_value('ID'));

						foreach my $primary_tag_key_level3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...

							if ( exists ($hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID} ) ){
								foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID}}) {
									embl_rendering($feature_level3);
								}
							}
						}
					}
				}
			}
		}
	}
}

sub embl_rendering {

	my ($feature)=@_;

	## check primary tag
	my $primary_tag = lc($feature->primary_tag);

	my @feature_list=["assembly_gap","C_region","CDS","centromere","D-loop","D_segment","exon","gap","gene","iDNA","intron","J_segment","LTR","mat_peptide","misc_binding","misc_difference","misc_feature","misc_recomb","misc_RNA","misc_structure","mobile_element","modified_base","mRNA","ncRNA","N_region","old_sequence","operon","oriT","polyA_site","precursor_RNA","prim_transcript","primer_bind","protein_bind","regulatory","repeat_region","rep_origin","rRNA","S_region","sig_peptide","source","stem_loop","STS","telomere","tmRNA","transit_peptide","tRNA","unsure","V_region","V_segment","variation","3'UTR","5'UTR"];

	foreach my $element (@feature_list){
		if(lc($element) =~ /$primary_tag/){
			$feature->$primary_tag = $element;
		}
		else{
			#repeat exception rule
			if( $primary_tag =~ /repeat/ ){
				$feature->$primary_tag = "repeat_region";
			}
			#utr5 exception rule
			elsif($primary_tag =~ /utr/ and ($primary_tag =~ /3/ or $primary_tag =~ /three/) ){
				$feature->$primary_tag = "3'UTR";
			}
			#utr5 exception rule
			elsif($primary_tag =~ /utr/ and ($primary_tag =~ /5/ or $primary_tag =~ /five/) ){
				$feature->$primary_tag = "5'UTR";
			}
			print "WARNING: this primary tag ".$primary_tag." is not recognized among those expected to be EMBL compliant. Please check it or create an exception rule.\n";
		}
	}
}

# check all the attribute to URI encode the values
sub uri_encode_omniscient {

	my ($hash_omniscient) = @_  ;

	#################
	# == LEVEL 1 == #
	#################
	foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
		foreach my $id_tag_key_level1 (keys %{$hash_omniscient->{'level1'}{$primary_tag_key_level1}}){
			uri_encode_one_feature($hash_omniscient->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1});

			#################
			# == LEVEL 2 == #
			#################
			foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

				if ( exists ($hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1} ) ){
					foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}) {
						uri_encode_one_feature($feature_level2);

						#################
						# == LEVEL 3 == #
						#################
						my $level2_ID  = lc($feature_level2->_tag_value('ID'));

						foreach my $primary_tag_key_level3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...

							if ( exists ($hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID} ) ){
								foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID}}) {
									uri_encode_one_feature($feature_level3);
								}
							}
						}
					}
				}
			}
		}
	}
}

sub uri_decode_omniscient {

	my ($hash_omniscient) = @_  ;

	#################
	# == LEVEL 1 == #
	#################
	foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
		foreach my $id_tag_key_level1 (keys %{$hash_omniscient->{'level1'}{$primary_tag_key_level1}}){
			uri_decode_one_feature($hash_omniscient->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1});

			#################
			# == LEVEL 2 == #
			#################
			foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

				if ( exists ($hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1} ) ){
					foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}) {
						uri_decode_one_feature($feature_level2);

						#################
						# == LEVEL 3 == #
						#################
						my $level2_ID  = lc($feature_level2->_tag_value('ID'));

						foreach my $primary_tag_key_level3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...

							if ( exists ($hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID} ) ){
								foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID}}) {
									uri_decode_one_feature($feature_level3);
								}
							}
						}
					}
				}
			}
		}
	}
}

# reencode in uri all value of all attributes of a feature
sub uri_encode_one_feature {

	my ($feature)=@_;

	my @list_tag = $feature->get_all_tags;

	foreach my $tag (@list_tag){
		my @values = $feature->get_tag_values($tag);
		$feature->remove_tag($tag);
		foreach my $val (@values){
			my $val_checked = uri_unescape($val);
			my $new_val = uri_escape($val_checked );
			$feature->add_tag_value($tag, $new_val);
		}
	}
}

sub uri_decode_one_feature {

	my ($feature)=@_;

	my @list_tag = $feature->get_all_tags;

	foreach my $tag (@list_tag){
		my @values = $feature->get_tag_values($tag);
		$feature->remove_tag($tag);
		foreach my $val (@values){
			my $new_val = uri_unescape($val);
			$feature->add_tag_value($tag, $new_val);
		}
	}
}


#				   +------------------------------------------------------+
#				   |+----------------------------------------------------+|
#				   || 			HANDLE OMNISCIENT => Fill / Modify		 ||
#				   |+----------------------------------------------------+|
#				   +------------------------------------------------------+



# omniscient is a hash containing a whole gXf file in memory sorted in a specific way (3 levels)
# If a feature/record already exists in omniscient_to_append, it will be replaced by the new one (If the new one content less features, the surnumerary ones are actually erased/removed)
sub fill_omniscient_from_other_omniscient_level1_id {

	my ($level_id_list, $hash_omniscient, $omniscient_to_append)=@_;

	#################
	# == LEVEL 1 == #
	#################
	foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
		foreach my $id_tag_key_level1 (@$level_id_list){ #select Id of the list
			if( exists_keys($hash_omniscient, ('level1', $primary_tag_key_level1, $id_tag_key_level1)) ){
				$omniscient_to_append->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1} = $hash_omniscient->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1}; # print feature

				#################
				# == LEVEL 2 == #
				#################
				foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
					if( exists_keys($hash_omniscient, ('level2', $primary_tag_key_level2, $id_tag_key_level1)) ){
						@{$omniscient_to_append->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}} = @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}};

						foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}) {

							#################
							# == LEVEL 3 == #
							#################
							my $level2_ID = lc($feature_level2->_tag_value('ID'));

							# remove feature from omniscient_to_append related to level2_ID. Like that those that are not anymore in hash_omniscient will not appear anymore at the end of the next step (#Now add the new L3 feature of level2_ID)
							foreach my $tag_l3 (keys %{$omniscient_to_append->{'level3'}}){ 
                                                                if( exists_keys($omniscient_to_append, ('level3', $tag_l3, $level2_ID)) ){
                                                                	delete $omniscient_to_append->{'level3'}{$tag_l3}{$level2_ID} ;
								}
                                                        }
							#Now add the new L3 feature of level2_ID
							foreach my $primary_tag_key_level3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
								if( exists_keys($hash_omniscient, ('level3', $primary_tag_key_level3, $level2_ID)) ){
									@{$omniscient_to_append->{'level3'}{$primary_tag_key_level3}{$level2_ID}} = @{$hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID}};
								}
							}
						}
					}
				}
			}
		}
	}
}

#append hash1 by hash2 accodingly with overlap parameter. Only non overlaping one will be kept
sub complement_omniscients_by_size {
	my ($omniscient1, $omniscient2, $size_min, $verbose)=@_;

	my %add_omniscient;

	$size_min=0 if (! $size_min);

	if(! $verbose){$verbose=0;}
	my $omniscient1_sorted = sort_by_seq($omniscient1);
	my $omniscient2_sorted = sort_by_seq($omniscient2);

	foreach my $locusID ( keys %{$omniscient2_sorted}){ # tag_l1 = gene or repeat etc...
		if ( exists_keys( $omniscient2_sorted, ( $locusID, 'level1') ) ){	
			foreach my $tag_l1 ( keys %{$omniscient2_sorted->{$locusID}{'level1'}} ) { 

				# Go through location from left to right ### !!
				foreach my $id1_l1 ( sort {$omniscient2_sorted->{$locusID}{'level1'}{$tag_l1}{$a}[1] <=> $omniscient2_sorted->{$locusID}{'level1'}{$tag_l1}{$b}[1] } keys %{$omniscient2_sorted->{$locusID}{'level1'}{$tag_l1}} ) {
					print "\nlets look at $id1_l1.\n" if ($verbose >= 3);
					my $take_it=1;
					my $location = $omniscient2_sorted->{$locusID}{'level1'}{$tag_l1}{$id1_l1};

					if( exists_keys($omniscient1_sorted, ($locusID,'level1',$tag_l1) ) ) {
						
						
						foreach my $id2_l1 ( sort {$omniscient1_sorted->{$locusID}{'level1'}{$tag_l1}{$a}[1] <=> $omniscient1_sorted->{$locusID}{'level1'}{$tag_l1}{$b}[1] } keys %{$omniscient1_sorted->{$locusID}{'level1'}{$tag_l1}} ) {
							
							my $location2 = $omniscient1_sorted->{$locusID}{'level1'}{$tag_l1}{$id2_l1}; # location hash2
														
							#If location_to_check start if over the end of the reference location, we stop
							if($location2->[1] > $location->[2]) {last;} 
							print "location_to_check start if over the end of the reference location.\n" if ($verbose >= 3);

							#If location_to_check end if inferior to the start of the reference location, we continue next
							if($location2->[2] < $location->[1]) {next;} 
							print "location_to_check start if inferior to the start of the reference location.\n" if ($verbose >= 3);

							# Let's check at Gene LEVEL
							if( location_overlap($location, $location2) ){ #location overlap at gene level check now level3
								#let's check at CDS level (/!\ id1_l1 is corresponding to id from $omniscient2)
								if(check_gene_overlap_at_CDSthenEXON($omniscient2, $omniscient1, $id1_l1, $id2_l1)){ #If contains CDS it has to overlap at CDS level, otherwise any type of feature level3 overlaping is sufficient to decide that they overlap
									print "$id2_l1 overlaps $id1_l1, we skip it.\n" if ($verbose >= 3);
									$take_it=undef; last;
								}
								print "$id2_l1 overlaps $id1_l1 overlap but not at CDS level.\n" if ($verbose >= 3);
							}
							else{
								print "$id2_l1 DO NOT OVERLAP $id1_l1.\n" if ($verbose >= 3);
							}
						}	
					}

					# We keep it because is not overlaping
					if($take_it){
						print "We take it : $id1_l1\n" if ($verbose >= 3);

						#look at size
						my $still_take_it=undef;
						foreach my $tag_l2 (keys %{$omniscient2->{'level2'}} ){
							if(exists_keys($omniscient2,('level2', $tag_l2, $id1_l1))){
								foreach my $feature_l2 ( @{$omniscient2->{'level2'}{$tag_l2}{$id1_l1}} ){
									my $id_l2 = $feature_l2->_tag_value('ID');
									foreach my $tag_l3 (keys %{$omniscient2->{'level3'}} ){
										if(exists_keys($omniscient2,('level3', 'cds', lc($id_l2)))){
											my $cds_size=0;
											foreach my $feature_l3 ( @{$omniscient2->{'level3'}{'cds'}{lc($id_l2)}} ){
												my $size=$feature_l3->end - $feature_l3->start +1;
												$cds_size += $size;
											}
											if($cds_size >= $size_min){
												$still_take_it=1;
												last;
											}
										}
									}
									last if $still_take_it;
								}
								last if $still_take_it;
							}
						}
						# We keep it because has size over threshold
						if ($still_take_it){
							#save level1
							$add_omniscient{'level1'}{$tag_l1}{$id1_l1} = $omniscient2->{'level1'}{$tag_l1}{$id1_l1};
							#save level2
							foreach my $tag_l2 (keys %{$omniscient2->{'level2'}} ){
								if(exists_keys($omniscient2,('level2', $tag_l2, $id1_l1))){
									# Add the level2 list data
									$add_omniscient{'level2'}{$tag_l2}{$id1_l1} = $omniscient2->{'level2'}{$tag_l2}{$id1_l1};
									# for each level2 get the level3 subfeatures 
									foreach my $feature_l2 ( @{$omniscient2->{'level2'}{$tag_l2}{$id1_l1}} ){
										my $id_l2 = $feature_l2->_tag_value('ID');
										#save level3
										foreach my $tag_l3 (keys %{$omniscient2->{'level3'}} ){
											if(exists_keys($omniscient2,('level3', $tag_l3, lc($id_l2)))){
												$add_omniscient{'level3'}{$tag_l3}{lc($id_l2)} = $omniscient2->{'level3'}{$tag_l3}{lc($id_l2)};
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	#Now populate hash1 with data from hash2
	merge_omniscients($omniscient1, \%add_omniscient);

	undef %add_omniscient;

	return $omniscient1;
}

#append hash1 by hash2 accodingly with overlap parameter. Only non overlaping one will be kept
sub complement_omniscients {
	my ($omniscient1, $omniscient2, $verbose)=@_;

	my %add_omniscient;
	
	if(! $verbose){$verbose=0;}
	my $omniscient1_sorted = sort_by_seq($omniscient1);
	my $omniscient2_sorted = sort_by_seq($omniscient2);

	foreach my $locusID ( keys %{$omniscient2_sorted}){ # tag_l1 = gene or repeat etc...
		if ( exists_keys( $omniscient2_sorted, ( $locusID, 'level1') ) ){	
			foreach my $tag_l1 ( keys %{$omniscient2_sorted->{$locusID}{'level1'}} ) { 

				# Go through location from left to right ### !!
				foreach my $id1_l1 ( sort {$omniscient2_sorted->{$locusID}{'level1'}{$tag_l1}{$a}[1] <=> $omniscient2_sorted->{$locusID}{'level1'}{$tag_l1}{$b}[1] } keys %{$omniscient2_sorted->{$locusID}{'level1'}{$tag_l1}} ) {
					print "\nlets look at $id1_l1.\n" if ($verbose >= 3);
					my $take_it=1;
					my $location = $omniscient2_sorted->{$locusID}{'level1'}{$tag_l1}{$id1_l1};

					if( exists_keys($omniscient1_sorted, ($locusID,'level1',$tag_l1) ) ) {
						
						
						foreach my $id2_l1 ( sort {$omniscient1_sorted->{$locusID}{'level1'}{$tag_l1}{$a}[1] <=> $omniscient1_sorted->{$locusID}{'level1'}{$tag_l1}{$b}[1] } keys %{$omniscient1_sorted->{$locusID}{'level1'}{$tag_l1}} ) {
							
							my $location2 = $omniscient1_sorted->{$locusID}{'level1'}{$tag_l1}{$id2_l1}; # location hash2
														
							#If location_to_check start if over the end of the reference location, we stop
							if($location2->[1] > $location->[2]) {last;} 
							print "location_to_check start if over the end of the reference location.\n" if ($verbose >= 3);

							#If location_to_check end if inferior to the start of the reference location, we continue next
							if($location2->[2] < $location->[1]) {next;} 
							print "location_to_check start if inferior to the start of the reference location.\n" if ($verbose >= 3);

							# Let's check at Gene LEVEL
							if( location_overlap($location, $location2) ){ #location overlap at gene level check now level3
								#let's check at CDS level (/!\ id1_l1 is corresponding to id from $omniscient2)
								if(check_gene_overlap_at_CDSthenEXON($omniscient2, $omniscient1, $id1_l1, $id2_l1)){ #If contains CDS it has to overlap at CDS level, otherwise any type of feature level3 overlaping is sufficient to decide that they overlap
									print "$id2_l1 overlaps $id1_l1, we skip it.\n" if ($verbose >= 3);
									$take_it=undef; last;
								}
								print "$id2_l1 overlaps $id1_l1 overlap but not at CDS level.\n" if ($verbose >= 3);
							}
							else{
								print "$id2_l1 DO NOT OVERLAP $id1_l1.\n" if ($verbose >= 3);
							}
						}	
					}

					if($take_it){
						print "We take it : $id1_l1\n" if ($verbose >= 3);
						#save level1
						$add_omniscient{'level1'}{$tag_l1}{$id1_l1} = $omniscient2->{'level1'}{$tag_l1}{$id1_l1};
						#save level2
						foreach my $tag_l2 (keys %{$omniscient2->{'level2'}} ){
							if(exists_keys($omniscient2,('level2', $tag_l2, $id1_l1))){
								# Add the level2 list data
								$add_omniscient{'level2'}{$tag_l2}{$id1_l1} = $omniscient2->{'level2'}{$tag_l2}{$id1_l1};
								# for each level2 get the level3 subfeatures 
								foreach my $feature_l2 ( @{$omniscient2->{'level2'}{$tag_l2}{$id1_l1}} ){
									my $id_l2 = $feature_l2->_tag_value('ID');
									#save level3
									foreach my $tag_l3 (keys %{$omniscient2->{'level3'}} ){
										if(exists_keys($omniscient2,('level3', $tag_l3, lc($id_l2)))){
											$add_omniscient{'level3'}{$tag_l3}{lc($id_l2)} = $omniscient2->{'level3'}{$tag_l3}{lc($id_l2)};
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	#Now populate hash1 with data from hash2
	merge_omniscients($omniscient1, \%add_omniscient);

	undef %add_omniscient;

	return $omniscient1;
}

# omniscient is a hash containing a whole gXf file in memory sorted in a specific way (3 levels)
# rename ID in hash_omniscient2 that already exist in hash_omniscient1
sub rename_ID_existing_in_omniscient {

	my ($hash_omniscient1, $hash_omniscient2, $verbose)=@_;

	if(! $verbose){$verbose=1;}

	my $hash_whole_IDs = get_all_IDs($hash_omniscient1);
	my $hash2_whole_IDs = get_all_IDs($hash_omniscient2);
	
	my %hash_miscCount;
	my $miscCount = \%hash_miscCount;
	my $resume_case=undef;

	#################
	# == LEVEL 1 == #
	#################
	foreach my $tag_l1 (keys %{$hash_omniscient2->{'level1'}}){ # tag_l1 = gene or repeat etc...
		foreach my $id_l1 (keys %{$hash_omniscient2->{'level1'}{$tag_l1}}){
			my $new_parent=undef;
			my $uID = $hash_omniscient2->{'level1'}{$tag_l1}{$id_l1}->_tag_value('ID');

			if ( exists ( $hash_whole_IDs->{$id_l1} ) ){
				$resume_case++;
				my $feature = $hash_omniscient2->{'level1'}{$tag_l1}{$id_l1};
				$uID = replace_by_uniq_ID( $feature, $hash_whole_IDs,  $hash2_whole_IDs, $miscCount);
				$hash_omniscient2->{'level1'}{$tag_l1}{lc($uID)} = delete $hash_omniscient2->{'level1'}{$tag_l1}{$id_l1}; # save feature
				$new_parent=1;
			}
			#################
			# == LEVEL 2 == #
			#################
			foreach my $tag_l2 (keys %{$hash_omniscient2->{'level2'}}){ # tag_l2 = mrna or mirna or ncrna or trna etc...
				
				if (exists_keys ($hash_omniscient2, ('level2', $tag_l2, $id_l1) ) ){ #Non present in hash2, we create a list with one element
					
					foreach my $feature_l2 ( @{$hash_omniscient2->{'level2'}{$tag_l2}{$id_l1}}) {
						
						my $new_parent_l2=undef;

						if($new_parent){
							create_or_replace_tag($feature_l2, 'Parent', $uID);
						}

						my $uID_l2 = $feature_l2->_tag_value('ID');
						my $id_l2 = lc($uID_l2);

						if ( exists ( $hash_whole_IDs->{$id_l2} ) ){

							$resume_case++;
							$uID_l2 = replace_by_uniq_ID($feature_l2, $hash_whole_IDs,  $hash2_whole_IDs, $miscCount);	
							$new_parent_l2=1;
						}
	
						#################
						# == LEVEL 3 == #
						#################
						foreach my $tag_l3 (keys %{$hash_omniscient2->{'level3'}}){ 

							if (exists_keys ($hash_omniscient2, ('level3', $tag_l3, $id_l2) ) ){ 
								
								foreach my $feature_l3 ( @{$hash_omniscient2->{'level3'}{$tag_l3}{$id_l2}}) {

									if($new_parent_l2){
										create_or_replace_tag($feature_l3, 'Parent', $uID_l2);
									}

									my $uID_l3 = $feature_l3->_tag_value('ID');
									my $id_l3 = lc($uID_l3);

									if ( exists ( $hash_whole_IDs->{$id_l2} ) ){
										$resume_case++;
										$uID_l3 = replace_by_uniq_ID($feature_l3, $hash_whole_IDs,  $hash2_whole_IDs, $miscCount);
										
									}
								}
								#save list feature level3
								if($new_parent_l2){
									$hash_omniscient2->{'level3'}{$tag_l3}{lc($uID_l2)} = delete $hash_omniscient2->{'level3'}{$tag_l3}{lc($id_l2)} ;
								}
							}
						}
					}
					#save list feature level2
					if($new_parent){
						$hash_omniscient2->{'level2'}{$tag_l2}{lc($uID)} = delete $hash_omniscient2->{'level2'}{$tag_l2}{lc($id_l1)};
					}
				}
			}
		}
	}
	print "we renamed $resume_case cases\n" if($verbose and $resume_case);

	return $hash_omniscient2;
}

# put data from hash_omniscient2 in hash_omniscient1
# Features are added even if they are identical. If they have similar name, new name will be given too.
sub merge_omniscients {
	# $hash_omniscient1 = omniscient to append !!!
	my ($hash_omniscient1, $hash_omniscient2, $hash_whole_IDs)=@_;
	
	if (! $hash_whole_IDs){
		$hash_whole_IDs = get_all_IDs($hash_omniscient1);
	}
	my $hash2_whole_IDs = get_all_IDs($hash_omniscient2);
	
	my %hash_miscCount;
	my $miscCount = \%hash_miscCount;

	#################
	# == LEVEL 1 == #
	#################
	foreach my $tag_l1 (keys %{$hash_omniscient2->{'level1'}}){ # tag_l1 = gene or repeat etc...
		foreach my $id_l1 (keys %{$hash_omniscient2->{'level1'}{$tag_l1}}){

			my $new_parent=undef;
			my $uID = $hash_omniscient2->{'level1'}{$tag_l1}{$id_l1}->_tag_value('ID');

			if ( ! exists ( $hash_whole_IDs->{$id_l1} ) ){
					$hash_omniscient1->{'level1'}{$tag_l1}{$id_l1} = $hash_omniscient2->{'level1'}{$tag_l1}{$id_l1}; # save feature level1
					$hash_whole_IDs->{$id_l1}++;
			}
			else{
				#print "INFO level1:  Parent $id_l1 already exist. We generate a new one to avoid collision !\n";
				my $feature = $hash_omniscient2->{'level1'}{$tag_l1}{$id_l1};
				$uID = replace_by_uniq_ID( $feature, $hash_whole_IDs,  $hash2_whole_IDs, $miscCount);
				$hash_omniscient1->{'level1'}{$tag_l1}{lc($uID)} = $hash_omniscient2->{'level1'}{$tag_l1}{$id_l1}; # save feature level1
				$new_parent=1;
			}

			#################
			# == LEVEL 2 == #
			#################
			foreach my $tag_l2 (keys %{$hash_omniscient2->{'level2'}}){ # tag_l2 = mrna or mirna or ncrna or trna etc...
				
				if (exists_keys ($hash_omniscient2, ('level2', $tag_l2, $id_l1) ) ){ #Non present in hash2, we create a list with one element
					
					foreach my $feature_l2 ( @{$hash_omniscient2->{'level2'}{$tag_l2}{$id_l1}}) {
						
						my $new_parent_l2=undef;
						if($new_parent){
							create_or_replace_tag($feature_l2, 'Parent', $hash_omniscient1->{'level1'}{$tag_l1}{lc($uID)}->_tag_value('ID'));
						}

						my $uID_l2 = $feature_l2->_tag_value('ID');
						my $id_l2 = lc($uID_l2);

						if ( exists ( $hash_whole_IDs->{$id_l2} ) ){
							
							#print "INFO level2:  Parent $id_l2 already exist. We generate a new one to avoid collision !\n";
							$uID_l2 = replace_by_uniq_ID($feature_l2, $hash_whole_IDs,  $hash2_whole_IDs, $miscCount);	
							$new_parent_l2=1;
						}
						else{$hash_whole_IDs->{$id_l2}++;}

						#################
						# == LEVEL 3 == #
						#################
						foreach my $tag_l3 (keys %{$hash_omniscient2->{'level3'}}){ 

							if (exists_keys ($hash_omniscient2, ('level3', $tag_l3, $id_l2) ) ){ 
								
								foreach my $feature_l3 ( @{$hash_omniscient2->{'level3'}{$tag_l3}{$id_l2}}) {

									if($new_parent_l2){
										create_or_replace_tag($feature_l3, 'Parent', $uID_l2);
									}

									my $uID_l3 = $feature_l3->_tag_value('ID');
									my $id_l3 = lc($uID_l3);

									if ( exists ( $hash_whole_IDs->{$id_l3} ) ){
									#	print "INFO level3:  Parent $id_l3 already exist. We generate a new one to avoid collision !\n";
										$uID_l3 = replace_by_uniq_ID($feature_l3, $hash_whole_IDs,  $hash2_whole_IDs, $miscCount);
									}
									else{$hash_whole_IDs->{$id_l3}++;}
								}
								#save list feature level3
								@{$hash_omniscient1->{'level3'}{$tag_l3}{lc($uID_l2)} } = @{ $hash_omniscient2->{'level3'}{$tag_l3}{$id_l2} };
							}
						}
					}
					#save list feature level2
					@{$hash_omniscient1->{'level2'}{$tag_l2}{lc($uID)}} = @{ $hash_omniscient2->{'level2'}{$tag_l2}{$id_l1} };
				}
			}
		}
	}
	return $hash_omniscient1, $hash_whole_IDs;
}

sub append_omniscient {

	my ($omniscient, $level1,$level2,$level3)=@_;

	foreach my $feature (@$level1){
		my $primaryTag = lc($feature->primary_tag);
		my $id  = lc($feature->_tag_value('ID'));

		if( ! exists_keys($omniscient, ('level1', $primaryTag, $id)) ){
			$omniscient->{"level1"}{$primaryTag}{$id}=$feature;
		}
	}
	foreach my $feature (@$level2){ # if exist, try to append the list
		my $primaryTag = lc($feature->primary_tag);
		my $parent_id  = lc($feature->_tag_value('Parent'));
		
		if( ! exists_keys($omniscient, ('level2', $primaryTag, $parent_id)) ){
			push(@{$omniscient->{"level2"}{$primaryTag}{$parent_id}}, $feature);###
		}
		else{ # append element in the list if not existing
			my $exist_in_list="no";
			my $id = lc($feature->_tag_value('ID'));
			
			foreach my $feature_original (@{$omniscient->{"level2"}{$primaryTag}{$parent_id}}){
				my $original_id = lc($feature_original->_tag_value('ID'));
				if ($original_id eq $id){
					$exist_in_list="yes"; last;
				}
			}
			if($exist_in_list eq "no"){ # feature doesnt exist in the feature list already present. So, we append it.
				push(@{$omniscient->{"level2"}{$primaryTag}{$parent_id}}, $feature)
			}
		}
	}
	foreach my $feature (@$level3){
		my $primaryTag = lc($feature->primary_tag);
		my $parent_id = lc($feature->_tag_value('Parent'));
		
		if( ! exists_keys($omniscient, ('level3', $primaryTag, $parent_id)) ){
			push(@{$omniscient->{"level3"}{$primaryTag}{$parent_id}}, $feature);
		}
		else{ # append element in the list if not existing
			my $exist_in_list="no";
			my $id = lc($feature->_tag_value('ID'));

			foreach my $feature_original (@{$omniscient->{"level3"}{$primaryTag}{$parent_id}}){

				my $original_id = lc($feature_original->get_tag_values('ID'));

				if ($original_id eq $id){
					$exist_in_list="yes"; last;
				}
			}
			if($exist_in_list eq "no"){ # feature doesnt exist in the feature list already present. So, we append it.
				push(@{$omniscient->{"level3"}{$primaryTag}{$parent_id}}, $feature)
			}
		}
	}
}

#				   +------------------------------------------------------+
#				   |+----------------------------------------------------+|
#				   || 			HANDLE OMNISCIENT => Remove				 ||
#				   |+----------------------------------------------------+|
#				   +------------------------------------------------------+

# omniscient is a hash containing a whole gXf file in memory sorted in a specific way (3 levels)
# Input: list of level1 id
#        omniscient
#
sub remove_omniscient_elements_from_level1_id_list {

	my ($hash_omniscient, $level_id_list) = @_  ;

	#################
	# == LEVEL 1 == #
	#################
	foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
		foreach my $id_tag_key_level1 (keys %{$hash_omniscient->{'level1'}{$primary_tag_key_level1}}){
			foreach my $level_id (@$level_id_list){
				if($id_tag_key_level1 eq lc($level_id)){

					#################
					# == LEVEL 2 == #
					#################
					foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

						if ( exists ($hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1} ) ){
							foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}) {

								#################
								# == LEVEL 3 == #
								#################
								my $level2_ID = lc($feature_level2->_tag_value('ID'));

								foreach my $primary_tag_key_level3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
									if ( exists ($hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID} ) ){
										delete $hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID} # delete level3
									}
								}
							}
							delete $hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1} # delete level2
						}
					}
				delete $hash_omniscient->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1}; # delete level1
				}
			}
		}
	}
}

# /!\XXX Has to be improved, we should loop over the feature list and extract the id_tag_key_level1 before to loop over the hash
# Input: list of level2 id
#        omniscient
#
sub remove_omniscient_elements_from_level2_feature_list {

	my ($hash_omniscient, $feature_list) = @_  ;

	#################
	# == LEVEL 2 == #
	#################
	foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
		foreach my $id_tag_key_level1 (keys %{$hash_omniscient->{'level2'}{$primary_tag_key_level2}}){
			if( exists_keys($hash_omniscient, ('level2', $primary_tag_key_level2, $id_tag_key_level1)) ){
				foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}) {
					my $level2_ID= lc($feature_level2->_tag_value('ID'));

					foreach my $feature (@$feature_list){
						my $feature_ID = lc($feature->_tag_value('ID'));
						my $feature_Parent_ID = lc($feature->_tag_value('Parent'));

						if($level2_ID eq $feature_ID){

							#################
							# == LEVEL 3 == #
							#################
							foreach my $primary_tag_key_level3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
								if( exists_keys($hash_omniscient, ('level3', $primary_tag_key_level3, $level2_ID)) ){
									delete $hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID} # delete level3
								}
							}
						my @id_concern_list=($feature_Parent_ID);
						my @id_list_to_remove=($feature_ID);
						my @list_tag_key=('all');
						remove_element_from_omniscient(\@id_concern_list, \@id_list_to_remove, $hash_omniscient, 'level2','false', \@list_tag_key);

							if( ! exists_keys($hash_omniscient, ('level2', $primary_tag_key_level2, $id_tag_key_level1)) ){
								#New new list was empty so l2 has been removed, we can now remove l1
								foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){
									if( exists_keys($hash_omniscient, ('level1', $primary_tag_key_level1, $feature_Parent_ID)) ){
										delete $hash_omniscient->{'level1'}{$primary_tag_key_level1}{$feature_Parent_ID}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

# After cleaning if nothing left attache to level 2 we removed it, and the same for level1
sub remove_omniscient_elements_from_level2_ID_list {

	my ($hash_omniscient, $ID_l2_list) = @_  ;

	#################
	# == LEVEL 2 == #
	#################
	foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
		foreach my $id_tag_key_level1 (keys %{$hash_omniscient->{'level2'}{$primary_tag_key_level2}}){
			if( exists_keys($hash_omniscient, ('level2', $primary_tag_key_level2, $id_tag_key_level1)) ){
				foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}) {
					my $level2_ID= lc($feature_level2->_tag_value('ID'));
					my $feature_Parent_ID = lc($feature_level2->_tag_value('Parent'));

					foreach my $feature_ID (@$ID_l2_list){
						$feature_ID = lc($feature_ID);
						
						if($level2_ID eq $feature_ID){

							#################
							# == LEVEL 3 == #
							#################
							foreach my $primary_tag_key_level3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
								if( exists_keys($hash_omniscient, ('level3', $primary_tag_key_level3, $level2_ID)) ){
									delete $hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID} # delete level3
								}
							}

							my @id_concern_list=($feature_Parent_ID);
							my @id_list_to_remove=($feature_ID);
							my @list_tag_key=('all');
							remove_element_from_omniscient(\@id_concern_list, \@id_list_to_remove, $hash_omniscient, 'level2','false', \@list_tag_key);

							if( ! exists_keys($hash_omniscient, ('level2', $primary_tag_key_level2, $id_tag_key_level1)) ){
								#New list was empty so l2 has been removed, we can now remove l1
								foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){
									if( exists_keys($hash_omniscient, ('level1', $primary_tag_key_level1, $feature_Parent_ID)) ){
										delete $hash_omniscient->{'level1'}{$primary_tag_key_level1}{$feature_Parent_ID}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

# remove value of hash from omniscient in level $level which have the tag incriminated
sub remove_tuple_from_omniscient {

	my ($id_list_to_remove, $hash_omniscient, $level, $bolean, $list_tag_key)=@_;

	# bolean true => we remove if in list_tag_key
	# bolean false => we remove if ti is not in list_tag_key
	my $remove;
	foreach my $tag_key  (keys %{$hash_omniscient->{$level}}){
		if($bolean eq 'true'){
			$remove="no";
		}else{$remove="yes";}
		foreach my $tag_key_to_match (@$list_tag_key){
			if ((lc($tag_key) eq  lc($tag_key_to_match)) and ($bolean eq 'true')){
				$remove="yes";
			}
			if ((lc($tag_key) eq  lc($tag_key_to_match)) and ($bolean eq 'false')){
				$remove="no";last;
			}
		}
		if ($remove eq 'yes'){
			foreach my $id_key  (keys %{$hash_omniscient->{$level}{$tag_key}}){
				foreach my $id_to_remove (@$id_list_to_remove){
					if(lc($id_to_remove) eq lc($id_key)){
						delete $hash_omniscient->{lc($level)}{lc($tag_key)}{lc($id_to_remove)}; #REMOVE THAT KEY-VALUE pair
					}
				}
			}
		}
	}
}

# from omniscient: remove feature from "feature list" of level2 or level3 with id present in $id_list_to_remove
# $id_concern = ID of parent we will check
sub remove_element_from_omniscient {

	my ($id_concern_list, $id_list_to_remove, $hash_omniscient, $level, $bolean, $list_tag_key)=@_;

	# bolean true => we remove if in list_tag_key
	# bolean false => we remove if is not in list_tag_key
	my $remove;
	#Check level and tag
	foreach my $tag_key  (keys %{$hash_omniscient->{$level}}){
		if($bolean eq 'true'){
			$remove="no";
		}else{$remove="yes";}
		foreach my $tag_key_to_match (@$list_tag_key){

			if ((lc($tag_key) eq  lc($tag_key_to_match)) and ($bolean eq 'true')){
				$remove="yes";
			}
			if ((lc($tag_key) eq  lc($tag_key_to_match)) and ($bolean eq 'false')){
				$remove="no";last;
			}
		}
		#Check feature id from list
		if ($remove eq 'yes'){
			foreach my $id_concern (@$id_concern_list){
				my $mustModifyList=undef;
				my @listok;

				if(exists_keys($hash_omniscient, ($level,$tag_key,lc($id_concern)))){
					foreach my $feature (@{$hash_omniscient->{$level}{$tag_key}{lc($id_concern)}}){
						my $id  = lc($feature->_tag_value('ID'));
						my $shouldremoveit=undef;

						foreach my $id_to_remove (@$id_list_to_remove){

							if(lc ($id_to_remove) eq $id){ # These feature is in list to remove
								$mustModifyList="yes"; $shouldremoveit="yes"; last;
							}
						}
						if(! $shouldremoveit){
							push(@listok, $feature);
						} # Feature not present in id_to_remove, we keep it in list.
					}
					if($mustModifyList){ # at least one feature has been removed from list. Save the new list
						if(@listok){
							@{$hash_omniscient->{$level}{$tag_key}{$id_concern}}=@listok;
						}
						else{ # The list is empty we could remove the key (otherwise we would have saved a emplty list)
							delete $hash_omniscient->{$level}{$tag_key}{$id_concern};
						}
					}
				}
			}
		}
	}
}

# $id_concern = ID of parent we will check
# Go trhough all the element (L1 or L2 list) and check if we find one with the specied tag attribute and value attribute. In that case we remove it from the list
sub remove_element_from_omniscient_attributeValueBased {

	my ($id_concern_list, $attributeValue, $attributeTag, $hash_omniscient, $level, $bolean, $list_tag_key)=@_;

	# bolean true => we remove if in list_tag_key
	# bolean false => we remove if is not in list_tag_key
	my $remove;
	#Check level and tag
	foreach my $tag_key  (keys %{$hash_omniscient->{$level}}){
		if($bolean eq 'true'){
			$remove="no";
		}else{$remove="yes";}
		foreach my $tag_key_to_match (@$list_tag_key){

			if ((lc($tag_key) eq  lc($tag_key_to_match)) and ($bolean eq 'true')){
				$remove="yes";
			}
			if ((lc($tag_key) eq  lc($tag_key_to_match)) and ($bolean eq 'false')){
				$remove="no";last;
			}
		}
		#Check feature id from list
		if ($remove eq 'yes'){
			foreach my $id_concern (@{$id_concern_list}){
				my $mustModifyList=undef;
				my @listok;

				if(exists_keys($hash_omniscient, ($level,$tag_key,lc($id_concern)))){
					foreach my $feature (@{$hash_omniscient->{$level}{$tag_key}{lc($id_concern)}}){
						my $id  = lc($feature->_tag_value('ID'));
						my $shouldremoveit=undef;

						if($feature->has_tag($attributeTag)){
							if( lc($feature->_tag_value($attributeTag)) eq lc($attributeValue) ){
								$mustModifyList="yes"; $shouldremoveit="yes";
							}
							else{push(@listok, $feature);}
						}
						else{push(@listok, $feature);}
					}
					if($mustModifyList){ # at least one feature has been removed from list. Save the new list
						@{$hash_omniscient->{$level}{$tag_key}{$id_concern}}=@listok;
					}
				}
			}
		}
	}
}

#				   +------------------------------------------------------+
#				   |+----------------------------------------------------+|
#				   || 			HANDLE OMNISCIENT => CREATE				 ||
#				   |+----------------------------------------------------+|
#				   +------------------------------------------------------+

sub create_omniscient {

	my ($level1,$level2,$level3)=@_;

	my $omniscient;

	foreach my $feature (@$level1){
		my $id = lc($feature->_tag_value('ID'));
		$omniscient->{"level1"}{lc($feature->primary_tag)}{$id}=$feature;
	}
	foreach my $feature (@$level2){
		my $id = lc($feature->_tag_value('Parent'));
		push(@{$omniscient->{"level2"}{lc($feature->primary_tag)}{$id}}, $feature);###
	}
	foreach my $feature (@$level3){
		my @parentList = lc( $feature->get_tag_values('Parent'));
		foreach my $id (@parentList){
			push(@{$omniscient->{"level3"}{lc($feature->primary_tag)}{$id}}, $feature);
		}
	}
	return $omniscient;
}

# This method will create a new omniscient from an omniscient of reference and a list of id level2
#$list_id_l2 has to be lower case
sub create_omniscient_from_idlevel2list{
	my ($omniscientref, $hash_mRNAGeneLink, $list_id_l2)=@_;

	my %omniscient_new;

	foreach my $id_l2 (@$list_id_l2){
		my  $id_l1 = $hash_mRNAGeneLink->{$id_l2};

		# ADD LEVEL1
		foreach my $tag_l1 (keys %{$omniscientref->{'level1'}}){
			if( exists ($omniscientref->{'level1'}{$tag_l1}{$id_l1}) ){
				$omniscient_new{'level1'}{$tag_l1}{$id_l1}=$omniscientref->{'level1'}{$tag_l1}{$id_l1};
				last;
			}
		}
		# ADD LEVEL2
		foreach my $tag_l2 (keys %{$omniscientref->{'level2'}}){
			if( exists ($omniscientref->{'level2'}{$tag_l2}{$id_l1}) ){
				foreach my $feature_l2 ( @{$omniscientref->{'level2'}{$tag_l2}{$id_l1}}){
					if(lc($feature_l2->_tag_value('ID')) eq $id_l2 ){
						push (@{$omniscient_new{'level2'}{$tag_l2}{$id_l1}}, $feature_l2);
						last;
					}
				}
			}
		}
		# ADD LEVEL3
		foreach my $tag_l3 (keys %{$omniscientref->{'level3'}}){
			if( exists ($omniscientref->{'level3'}{$tag_l3}{$id_l2}) ){
				foreach my $feature_l3 ( @{$omniscientref->{'level3'}{$tag_l3}{$id_l2}}){
					push (@{$omniscient_new{'level3'}{$tag_l3}{$id_l2}}, $feature_l3);
				}
			}
		}
	}
	return \%omniscient_new;
}

# @Purpose: filter an omniscient to return a new omnicient containing only data related by the list of level1 IDs
# @input: 1 =>  omniscient hash reference
# @output 1 =>  omniscient hash reference
sub subsample_omniscient_from_level1_id_list {

	my ($hash_omniscient, $level_id_list) = @_  ;

	my %new_hash;

	#################
	# == LEVEL 1 == #
	#################
	foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...

		foreach my $id_tag_key_level1_raw (@$level_id_list){
			my $id_tag_key_level1 = lc($id_tag_key_level1_raw);
			if(exists ($hash_omniscient->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1})){

				$new_hash{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1} = delete $hash_omniscient->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1};

				#################
				# == LEVEL 2 == #
				#################
				foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
					if ( exists ($hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1} ) ){
						foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}) {
							my $level2_ID = lc($feature_level2->_tag_value('ID'));

							#################
							# == LEVEL 3 == #
							#################
							foreach my $primary_tag_key_level3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
								if ( exists ($hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID} ) ){
									$new_hash{'level3'}{$primary_tag_key_level3}{$level2_ID} = delete $hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID};
								}
							}
						}
						$new_hash{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1} = delete $hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1};
					}
				}
			}
		}
	}
	return \%new_hash;
}

#				   +------------------------------------------------------+
#				   |+----------------------------------------------------+|
#				   || 					Miscenaleous					 ||
#				   |+----------------------------------------------------+|
#				   +------------------------------------------------------+


# INPUT: feature object, String tag, String or Array ref;
# Output: None
sub create_or_replace_tag{

	my ($feature, $tag, $value)=@_;

	if ($feature->has_tag($tag) ) {
			$feature->remove_tag($tag);
			if(ref($value) eq "ARRAY"){
				$feature->add_tag_value($tag,@{$value});
			}
			else{
        		$feature->add_tag_value($tag,$value);
        	}
	}
	else{
		if(ref($value) eq "ARRAY"){
			$feature->add_tag_value($tag,@{$value});
		}
		else{
			$feature->add_tag_value($tag,$value);
		}
	}
}

# frame explanation
#0 indicates that the feature begins with a whole codon at the 5' most base. 1 means that there is one extra base (the third base of a codon) before the first whole codon
#and 2 means that there are two extra bases (the second and third bases of the codon) before the first codon.
sub fil_cds_frame {

	my ($hash_omniscient)=@_;

	foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

		foreach my $id_tag_key_level1 (keys %{$hash_omniscient->{'level2'}{$primary_tag_key_level2}}) {

			foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}) {

				my @temp = $feature_level2->get_tag_values('ID');
				my $level2_ID = lc(shift @temp);

				# == LEVEL 3 == #
				if ( exists_keys($hash_omniscient,('level3','cds',$level2_ID) ) ){
					my $strand=$feature_level2->strand;
					my @cds_list;
					if(($feature_level2->strand eq "+") or ($feature_level2->strand eq "1")){
						@cds_list=sort {$a->start <=> $b->start}  @{$hash_omniscient->{'level3'}{'cds'}{$level2_ID}};
					}else{
						@cds_list=sort {$b->start <=> $a->start}  @{$hash_omniscient->{'level3'}{'cds'}{$level2_ID}};
					}


					my $extra_base_start_cds=0;
					my $frame=0;
					foreach my $cds_feature ( @cds_list) {
						$cds_feature->frame($frame);
						my $cds_length=$cds_feature->end-$cds_feature->start +1;
						$frame=(3-(($cds_length-$frame)%3))%3; #second modulo allows to avoid the frame with 3. Instead we have 0.
					}
				}
			}
		}
	}
}

sub info_omniscient {
	my ($hash_omniscient)=@_;

	my %resu;

	foreach my $tag (keys %{$hash_omniscient->{'level1'}}){
    	my $nb=keys %{$hash_omniscient->{'level1'}{$tag}};
    	$resu{$tag}=$nb;
	}

	foreach my $level (keys %{$hash_omniscient}){
  		if ($level ne 'level1'){
    			foreach my $tag (keys %{$hash_omniscient->{$level}}){
      				foreach my $id (keys %{$hash_omniscient->{$level}{$tag}}){
        				my $nb=$#{$hash_omniscient->{$level}{$tag}{$id}}+1;
						if(exists_keys(\%resu,($tag))){
						        $resu{$tag}=$resu{$tag}+$nb;
						}
						else{
							$resu{$tag}=$nb;
						}
     				}
    			}
  		}
	}
	foreach my $tag (keys %resu){
		print "There is $resu{$tag} $tag\n";
	}
}

#check if reference exists in hash. Deep infinite : hash{a} or hash{a}{b} or hash{a}{b}{c}, etc.
# usage example: exists_keys($hash_omniscient,('level3','cds',$level2_ID)
sub exists_keys {
    my ($hash, @keys) = @_;

    for my $key (@keys) {
    	if (ref $hash ne 'HASH' or ! exists $hash->{$key}) {
    		return '';
    	}
		$hash = $hash->{$key};
    }
    return 1;
}

# omniscient is a hash containing a whole gXf file in memory sorted in a specific way (3 levels)
# This method group all features of a seq_id together.
sub group_features_from_omniscient {

	my ($hash_omniscient) = @_  ;

	my %group;
	#################
	# == LEVEL 1 == #
	#################
	foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
		my $key;
		foreach my $id_tag_key_level1 (keys %{$hash_omniscient->{'level1'}{$primary_tag_key_level1}}){
			my $feature_l1=$hash_omniscient->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1};
			my $seq_id=$feature_l1->seq_id;
			$key="$primary_tag_key_level1$id_tag_key_level1";
			push(@{$group{$seq_id}{$key}}, $hash_omniscient->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1});
			#################
			# == LEVEL 2 == #
			#################
			foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

				if ( exists ($hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1} ) ){
					foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}) {
						push(@{$group{$seq_id}{$key}}, $feature_level2);
						#################
						# == LEVEL 3 == #
						#################
						my $level2_ID = lc($feature_level2->_tag_value('ID'));

						############
						# THEN ALL THE REST
						foreach my $primary_tag_key_level3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...

							if ( exists ($hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID} ) ){
								foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID}}) {
									push(@{$group{$seq_id}{$key}}, $feature_level3);
								}
							}
						}
					}
				}
			}
		}
	}
	return \%group;
}

# omniscient is a hash containing a whole gXf file in memory sorted in a specific way (3 levels)
# This method group all level1 features of a seq_id together.
# hash{seq_id} = @(feature1, feature2 ...)
sub group_l1features_from_omniscient {

	my ($hash_omniscient) = @_  ;

	my %group;
	#################
	# == LEVEL 1 == #
	#################
	foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
		foreach my $id_tag_key_level1 (keys %{$hash_omniscient->{'level1'}{$primary_tag_key_level1}}){
			my $feature_l1=$hash_omniscient->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1};
			my $seq_id=$feature_l1->seq_id;
			push(@{$group{$seq_id}}, $feature_l1);

		}
	}
	return \%group;
}

# omniscient is a hash containing a whole gXf file in memory sorted in a specific way (3 levels)
# This method group all level1 features of a seq_id together.
# hash{seq_id} = @(id1, id2 ...)
sub group_l1IDs_from_omniscient {

	my ($hash_omniscient) = @_  ;

	my %group;
	#################
	# == LEVEL 1 == #
	#################
	foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
		foreach my $id_tag_key_level1 (keys %{$hash_omniscient->{'level1'}{$primary_tag_key_level1}}){
			my $feature_l1=$hash_omniscient->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1};
			my $seq_id=$feature_l1->seq_id;
			push(@{$group{$seq_id}}, lc($feature_l1->_tag_value('ID')));

		}
	}
	return \%group;
}

sub get_feature_l2_from_id_l2_l1 {
	my ($hash_omniscient, $id_l2, $id_l1) = @_  ;
	foreach my $tag_l2 (keys %{$hash_omniscient->{'level2'}}){
		if(exists ($hash_omniscient->{'level2'}{$tag_l2}{$id_l1})){
			foreach my $feature (@{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}}){
				if ( lc($feature->_tag_value('ID')) eq lc($id_l2) ) {
					return $feature
				}
			}
		}
		else{print "element level2 $tag_l2 $id_l1 doesnt exists in omniscient\n";}
	}
}

#extract sequences form list of cds features in a fasta db
# return a Bio::Seq object
sub extract_cds_sequence {
	my ($feature_list, $db)=@_;

	my $sequence="";
	foreach my $feature (sort {$a->start <=> $b->start} @$feature_list){
		$sequence .= $db->subseq($feature->seq_id,$feature->start,$feature->end);
	}
	my $seq  = Bio::Seq->new( '-format' => 'fasta' , -seq => $sequence);
	if($feature_list->[0]->strand eq "-1" or $feature_list->[0]->strand eq "-"){
		$seq=$seq->revcom;
	}
	return $seq ;
}

# @Purpose: from a omniscient and a gene_id, will get back the extrem value for start and end
# @input: 2 => hash(omniscient), string(gene identifier)
# @output: 2 => integer(extrem start position), integer(extrem end position)
sub get_longest_cds_start_end {
  my  ($hash_omniscient,$gene_id)=@_;
  my $resu_start=100000000000;
  my $resu_end=0;

  #check full CDS for each mRNA
  foreach my $mrna_feature (@{$hash_omniscient->{'level2'}{'mrna'}{lc($gene_id)}}){
    my $mrna_id = lc($mrna_feature->_tag_value('ID'));
    my $extrem_start=100000000000;
    my $extrem_end=0;

    #check all cds pieces
    foreach my $cds_feature (@{$hash_omniscient->{'level3'}{'cds'}{$mrna_id}}){
      if ($cds_feature->start < $extrem_start){
        $extrem_start=$cds_feature->start;
      }
      if($cds_feature->end > $extrem_end){
              $extrem_end=$cds_feature->end ;
      }
    }

    if($extrem_start < $resu_start){
        $resu_start=$extrem_start;
    }
    if($extrem_end > $resu_end){
      $resu_end=$extrem_end;
    }
  }
  return $resu_start,$resu_end;
}

# @Purpose: Filter the mRNA to keep only the one containing the longest CDS per gene
# @input: 1 => hash(omniscient hash)
# @output: list of id level2
sub get_longest_cds_level2{
  my ($hash_omniscient)= @_;

  my @list_id_l2;

  #################
  # == LEVEL 1 == #
  #################
  foreach my $primary_tag_l1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
    foreach my $id_tag_l1 (keys %{$hash_omniscient->{'level1'}{$primary_tag_l1}}){

      #################
      # == LEVEL 2 == #
      #################
      foreach my $primary_tag_l2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
        if ( exists ($hash_omniscient->{'level2'}{$primary_tag_l2}{$id_tag_l1} ) ){

          #check if there is isoforms
          ###########################

          #take only le longest
          if ($#{$hash_omniscient->{'level2'}{$primary_tag_l2}{$id_tag_l1}} > 0){
            my $longestL2 ="";
            my $longestCDSsize = 0;
            my $longestEXONsize = 0;

            foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_l2}{$id_tag_l1}}) {

	            my $level2_ID =   lc($feature_level2->_tag_value('ID') ) ;
	            if ( exists_keys( $hash_omniscient, ('level3','cds',$level2_ID ) ) ) {
	                my $cdsSize=0;
	                foreach my $cds ( @{$hash_omniscient->{'level3'}{'cds'}{$level2_ID}} ) { # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
	                  $cdsSize += ( $cds->end - $cds->start + 1 );
	                }

	                if($cdsSize > $longestCDSsize ){
	                  $longestL2 = $level2_ID;
	                  $longestCDSsize = $cdsSize;
	                }
	            }
	            elsif ( exists_keys( $hash_omniscient, ('level3','exon',$level2_ID ) ) ) {
					my $exonSize=0;
	                foreach my $exon ( @{$hash_omniscient->{'level3'}{'exon'}{$level2_ID}} ) { # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
	                  $exonSize += ( $exon->end - $exon->start + 1 );
	                }

	                if($exonSize > $longestEXONsize ){
	                  $longestL2 = $level2_ID;
	                  $longestEXONsize = $exonSize;
	                }
	            }
	            else{
	            	warn "WARNING get_longest_cds_level2: NO exon or cds to select the longest l2 for $id_tag_l1 l1 ! We will take one randomly ! @\n";
	            	$longestL2 = $level2_ID;
	        	}
            }
            push @list_id_l2,$longestL2; # push id of the longest
          }
          else{ #take it only of cds exits
            my $level2_ID =  lc(@{$hash_omniscient->{'level2'}{$primary_tag_l2}{$id_tag_l1}}[0]->_tag_value('ID')) ;
            if (exists_keys( $hash_omniscient, ('level3','cds', $level2_ID ) ) ){
              push @list_id_l2, $level2_ID; # push the only one existing
            } 
          }
        }
      }
    }
  }

  return \@list_id_l2;
}

# @Purpose: Counter the number of feature level in an omniscient
# @input: 1 => hash(omniscient hash)
# @output: integer
sub nb_feature_level1 {

  my ($omniscient)=@_;
  my $resu=0;
	foreach my $tag_level1 (keys %{$omniscient->{'level1'}}){
		$resu += (keys %{$omniscient->{'level1'}{$tag_level1}})
	}
	return $resu;
}

# @Purpose: get all the ID present in an omniscient
# @input: 1 => hash(omniscient hash)
# @output: hash of the whole IDs
sub get_all_IDs{
	my ($omniscient)=@_;

	my %whole_IDs;

	#################
	# == LEVEL 1 == #
	#################
	foreach my $primary_tag_l1 (keys %{$omniscient->{'level1'}}){ # primary_tag_l1 = gene or repeat etc...
		foreach my $id_l1 (keys %{$omniscient->{'level1'}{$primary_tag_l1}}){
			$whole_IDs{$id_l1}++;
		}
	}
	#################
	# == LEVEL 2 == #
	#################
	foreach my $primary_tag_l2 (keys %{$omniscient->{'level2'}}){ # primary_tag_l2 = mrna or mirna or ncrna or trna etc...
		foreach my $id_l1 ( keys %{$omniscient->{'level2'}{$primary_tag_l2}}) {
			foreach my $feature_level2 ( @{$omniscient->{'level2'}{$primary_tag_l2}{$id_l1}}) {
				my $level2_ID  = lc($feature_level2->_tag_value('ID'));		
				$whole_IDs{$level2_ID}++;
			}
		}
	}
	#################
	# == LEVEL 3 == #
	#################
	foreach my $primary_tag_l3 (keys %{$omniscient->{'level3'}}){ # primary_tag_l3 = cds or exon or start_codon or utr etc...
		foreach my $level2_ID ( keys %{$omniscient->{'level3'}{$primary_tag_l3}}) {
			foreach my $feature_level3 ( @{$omniscient->{'level3'}{$primary_tag_l3}{$level2_ID}}) {
				my $level3_ID  = lc($feature_level3->_tag_value('ID'));		
				$whole_IDs{$level3_ID}++;
			}
		}
	}
	return \%whole_IDs;
}

# @Purpose: Replace ID by Uniq ID and modify all parent attribute of child feature to stay in line with the modification
# @input: 4 => feature objetc, hash of ids, hash of ids, hash of feature counted to give more rapidly a name 
# @output: uniq ID
sub replace_by_uniq_ID{
	my ($feature, $hash_whole_IDs, $hash2_whole_IDs, $miscCount) = @_;

	my $id = $feature->_tag_value('ID');
	my $prefix = "IDmodified";
	my $key;

	if($prefix){
		$key=$prefix."-".lc($feature->primary_tag);
	}
	else{
		$key=lc($feature->primary_tag);
	}

	my $uID=$id;
	while( exists_keys($hash_whole_IDs, (lc($uID)) ) or exists_keys($hash2_whole_IDs, (lc($uID)) ) ){	 #loop until we found an uniq tag	
		$miscCount->{$key}++;
		$uID = $key."-".$miscCount->{$key};
	}

	#push the new ID	
	$hash_whole_IDs->{lc($uID)}=$id;

	# modify the feature ID with the correct one chosen
	create_or_replace_tag($feature,'ID', $uID); #modify ID to replace by parent value

	#Now repercute this modification to the subfeatures
	return $uID;
}

#				   +------------------------------------------------------+
#				   |+----------------------------------------------------+|
#				   || 			MANIPULATION AT OMNISCIENT LEVEL1/2/3	 ||
#				   |+----------------------------------------------------+|
#				   +------------------------------------------------------+

# @Purpose: Check 2 lists of feature L2 and remove the identical ones from the second list.
# @input: 4 =>  omniscient Hash reference, list1 reference of L2 features,  list2 reference of L2 features, verbose option for debug
# @output: list2 minus all the feature identical to one of the list1 feature
sub keep_only_uniq_from_list2{
	my ($omniscient, $list1_l2, $list2_l2, $verbose)= @_;

	my @new_list2;
	my $keep = 1;

	foreach my $feature2 ( @{$list2_l2} ){
		foreach my $feature1 ( @{$list1_l2} ){	
			if(l2_identical($omniscient, $feature1, $feature2, $verbose )){
				$keep = undef; last;
			}
		}
		if($keep){
			push(@new_list2, $feature2);
		}
		else{ # We dont keep the l2 feature so we have to remove all related features
			remove_l2_related_feature($omniscient, $feature2, $verbose);
		}
	}
	return \@new_list2;
}

# check if l2 are identical. So look recursively at the level under.
# return 1 if identical
sub l2_identical{
	my ($omniscient, $feature1_l2, $feature2_l2, $verbose)= @_;
	my $result=1;

	my $id1_l2 = lc($feature1_l2->_tag_value('ID') );
	my $id2_l2 = lc($feature2_l2->_tag_value('ID') );

	foreach my $l3_type (keys %{$omniscient->{'level3'}} ){
		if(exists_keys($omniscient,('level3', $l3_type, $id1_l2))){
			if(exists_keys($omniscient,('level3', $l3_type, $id2_l2))){

				if(scalar @{$omniscient->{'level3'}{$l3_type}{$id1_l2}} ==  scalar @{$omniscient->{'level3'}{$l3_type}{$id2_l2}}){

					foreach my $feature1_level3 ( sort {$a->start <=> $b->start} @{$omniscient->{'level3'}{$l3_type}{$id1_l2}}) {
						
						my $identik = undef;
						foreach my $feature2_level3 ( sort {$a->start <=> $b->start} @{$omniscient->{'level3'}{$l3_type}{$id2_l2}}) {
						
							if( ($feature1_level3->start == $feature2_level3->start) and ($feature1_level3->end == $feature2_level3->end) ){
								$identik=1;	
							}
						}
						if(! $identik){
							return undef;
						}
					}
				}
				else{return undef;} # Not same number of features. Cannot be identical
			}
			else{return undef;}
		}
		else{ 
			if(exists_keys($omniscient,('level3', $l3_type, $id2_l2))){ # $id1_l2 do not have l3 but $id2_l2 has !
				return undef;
			}
		}
	}
	print "The isoforms $id1_l2 and $id2_l2 are identical\n" if ($verbose >= 2 and $result);
    return $result;
}

#
#
#Remove everything related to l2. But not itself ... why ?
sub remove_l2_related_feature{
	my ($omniscient, $feature2, $verbose) = @_;

	my $l1_id = lc($feature2->_tag_value('Parent'));
	my $l2_id = lc($feature2->_tag_value('ID'));

	#remove level 1 feature
	foreach my $tag (keys %{$omniscient->{'level1'}}){
		if(exists_keys($omniscient, ('level1', $tag, $l1_id))){
			delete $omniscient->{'level1'}{$tag}{$l1_id};
			last;
		}
	}
	foreach my $tag (keys %{$omniscient->{'level3'}}){
		if(exists_keys($omniscient, ('level3', $tag, $l2_id))){
			delete $omniscient->{'level3'}{$tag}{$l2_id};
		}
	}
}


# @Purpose: Convert the omniscient to be constistent with the ensembl gff format
# @input: 1 =>  omniscient Hash reference
# @output 1 =>  omniscient Hash reference
sub convert_omniscient_to_ensembl_style{
	my ($omniscient) = @_;

	convert_3th_column($omniscient, "mRNA", "transcript");
	add_exon_id($omniscient, "ID");
	add_transcript_id($omniscient, "ID");
	add_gene_id($omniscient, "ID");
}

# @Purpose: Convert the 3th column (feature type) from old value to new value
# @input: 3 =>  omniscient Hash reference, String = featurey type original, String = new feature type
# @output none =>  The hash itself is modified
sub convert_3th_column{
	my ($omniscient, $original, $new) = @_;

	foreach my $primary_tag_l1 (keys %{$omniscient->{'level1'}}){ # primary_tag_l1 = gene or repeat etc...
		foreach my $id_l1 (keys %{$omniscient->{'level1'}{$primary_tag_l1}}){
			
			if( lc($omniscient->{'level1'}{$primary_tag_l1}{$id_l1}->primary_tag) eq lc($original) ){
				$omniscient->{'level1'}{$primary_tag_l1}{$id_l1}->primary_tag($new);
			}

			#################
			# == LEVEL 2 == #
			#################
			foreach my $primary_tag_l2 (keys %{$omniscient->{'level2'}}){ # primary_tag_l2 = mrna or mirna or ncrna or trna etc...
				if ( exists ($omniscient->{'level2'}{$primary_tag_l2}{$id_l1} ) ){
					foreach my $feature_level2 ( @{$omniscient->{'level2'}{$primary_tag_l2}{$id_l1}}) {

						if(lc($feature_level2->primary_tag) eq lc($original) ){
							$feature_level2->primary_tag($new);
						}
						#################
						# == LEVEL 3 == #
						#################
						my $level2_ID  = lc($feature_level2->_tag_value('ID'));

						foreach my $primary_tag_l3 (keys %{$omniscient->{'level3'}}){ # primary_tag_l3 = cds or exon or start_codon or utr etc...
							if ( exists ($omniscient->{'level3'}{$primary_tag_l3}{$level2_ID} ) ){
								foreach my $feature_level3 ( @{$omniscient->{'level3'}{$primary_tag_l3}{$level2_ID}}) {
									if(lc($feature_level3->primary_tag) eq lc($original) ){
										$feature_level3->primary_tag($new);
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

# @Purpose: Copy past an attribute and change its tag
# @input: 3 =>  omniscient Hash reference, String = attribute tag original, String = attribute tag new
# @output none =>  The hash itself is modified
sub create_attribute_from_existing_attribute{
	my ($omniscient, $original_attribute, $new_attribute) = @_;

	foreach my $primary_tag_l1 (keys %{$omniscient->{'level1'}}){ # primary_tag_l1 = gene or repeat etc...
		foreach my $id_l1 (keys %{$omniscient->{'level1'}{$primary_tag_l1}}){
			my $feature_l1 = $omniscient->{'level1'}{$primary_tag_l1}{$id_l1};

			if( $feature_l1->has_tag($original_attribute)  and ! $feature_l1->has_tag($new_attribute)  ) {
				create_or_replace_tag($feature_l1,$new_attribute, $feature_l1->get_tag_values($original_attribute));
			}

			#################
			# == LEVEL 2 == #
			#################
			foreach my $primary_tag_l2 (keys %{$omniscient->{'level2'}}){ # primary_tag_l2 = mrna or mirna or ncrna or trna etc...
				if ( exists ($omniscient->{'level2'}{$primary_tag_l2}{$id_l1} ) ){
					foreach my $feature_l2 ( @{$omniscient->{'level2'}{$primary_tag_l2}{$id_l1}}) {

						if( $feature_l2->has_tag($original_attribute)  and ! $feature_l2->has_tag($new_attribute) ){
							create_or_replace_tag($feature_l2,$new_attribute, $feature_l2->get_tag_values($original_attribute));
						}
						#################
						# == LEVEL 3 == #
						#################
						my $level2_ID  = lc($feature_l2->_tag_value('ID'));

						foreach my $primary_tag_l3 (keys %{$omniscient->{'level3'}}){ # primary_tag_l3 = cds or exon or start_codon or utr etc...
							if ( exists ($omniscient->{'level3'}{$primary_tag_l3}{$level2_ID} ) ){
								foreach my $feature_l3 ( @{$omniscient->{'level3'}{$primary_tag_l3}{$level2_ID}}) {
									if( $feature_l3->has_tag($original_attribute) and ! $feature_l3->has_tag($new_attribute) ){
										create_or_replace_tag($feature_l3, $new_attribute, $feature_l3->get_tag_values($original_attribute));
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

# @Purpose: Create a locus tag for all feature L1 and share it with all children features
# @input: 3 =>  omniscient Hash reference, String = attribute tag original to  use as locus tag, String = locus_tag attribute tag
# @output none =>  The hash itself is modified
sub create_locus_tag{
	my ($omniscient, $original_attribute, $locus_tag) = @_;

	foreach my $primary_tag_l1 (keys %{$omniscient->{'level1'}}){ # primary_tag_l1 = gene or repeat etc...
		foreach my $id_l1 (keys %{$omniscient->{'level1'}{$primary_tag_l1}}){
			my $feature_l1 = $omniscient->{'level1'}{$primary_tag_l1}{$id_l1};

			my $locus_tag_value=undef;
			if( $feature_l1->has_tag($original_attribute)  and ! $feature_l1->has_tag($locus_tag)  ) {
				$locus_tag_value =  $feature_l1->_tag_value($original_attribute);
				create_or_replace_tag($feature_l1, $locus_tag, $locus_tag_value);
			}
			else{
				$locus_tag_value =  $feature_l1->_tag_value($locus_tag);
			}

			#################
			# == LEVEL 2 == #
			#################
			foreach my $primary_tag_l2 (keys %{$omniscient->{'level2'}}){ # primary_tag_l2 = mrna or mirna or ncrna or trna etc...
				if ( exists ($omniscient->{'level2'}{$primary_tag_l2}{$id_l1} ) ){
					foreach my $feature_l2 ( @{$omniscient->{'level2'}{$primary_tag_l2}{$id_l1}}) {

						if( $feature_l2->has_tag($original_attribute)  and ! $feature_l2->has_tag($locus_tag) ){
							create_or_replace_tag($feature_l2,$locus_tag, $locus_tag_value);
						}
						#################
						# == LEVEL 3 == #
						#################
						my $level2_ID  = lc($feature_l2->_tag_value('ID'));

						foreach my $primary_tag_l3 (keys %{$omniscient->{'level3'}}){ # primary_tag_l3 = cds or exon or start_codon or utr etc...
							if ( exists ($omniscient->{'level3'}{$primary_tag_l3}{$level2_ID} ) ){
								foreach my $feature_l3 ( @{$omniscient->{'level3'}{$primary_tag_l3}{$level2_ID}}) {
									if( $feature_l3->has_tag($original_attribute) and ! $feature_l3->has_tag($locus_tag) ){
										create_or_replace_tag($feature_l3, $locus_tag, $locus_tag_value);
									}
								}
							}
						}
					}
				}
			}
		}
	}
}


# @Purpose: add exon_id to exon feature if does not exist.
# @input: 3 =>  omniscient Hash reference, String = attribute to use to create exon_id
# @output none =>  The hash itself is modified
sub add_exon_id{
	my ($omniscient, $original_attribute) = @_;

	foreach my $primary_tag_l1 (keys %{$omniscient->{'level1'}}){ # primary_tag_l1 = gene or repeat etc...
		foreach my $id_l1 (keys %{$omniscient->{'level1'}{$primary_tag_l1}}){
			#################
			# == LEVEL 2 == #
			#################
			foreach my $primary_tag_l2 (keys %{$omniscient->{'level2'}}){ # primary_tag_l2 = mrna or mirna or ncrna or trna etc...
				if ( exists ($omniscient->{'level2'}{$primary_tag_l2}{$id_l1} ) ){
					foreach my $feature_l2 ( @{$omniscient->{'level2'}{$primary_tag_l2}{$id_l1}}) {
						#################
						# == LEVEL 3 == #
						#################
						my $level2_ID  = lc($feature_l2->_tag_value('ID'));
						foreach my $primary_tag_l3 (keys %{$omniscient->{'level3'}}){ # primary_tag_l3 = cds or exon or start_codon or utr etc...
							if ( exists ($omniscient->{'level3'}{$primary_tag_l3}{$level2_ID} ) ){
								foreach my $feature_l3 ( @{$omniscient->{'level3'}{$primary_tag_l3}{$level2_ID}}) {
									
									if( $feature_l3->has_tag($original_attribute) and ! $feature_l3->has_tag('exon_id') ){
										create_or_replace_tag($feature_l3, "exon_id", $feature_l3->_tag_value($original_attribute) );
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

# @Purpose: add transcript_id to all feature l2 if does not exist.
# @input: 3 =>  omniscient Hash reference, String = attribute to use to create transcript_id
# @output none =>  The hash itself is modified
sub add_transcript_id{
	my ($omniscient, $original_attribute) = @_;

	foreach my $primary_tag_l1 (keys %{$omniscient->{'level1'}}){ # primary_tag_l1 = gene or repeat etc...
		foreach my $id_l1 (keys %{$omniscient->{'level1'}{$primary_tag_l1}}){
			#################
			# == LEVEL 2 == #
			#################
			foreach my $primary_tag_l2 (keys %{$omniscient->{'level2'}}){ # primary_tag_l2 = mrna or mirna or ncrna or trna etc...
				if ( exists ($omniscient->{'level2'}{$primary_tag_l2}{$id_l1} ) ){
					foreach my $feature_l2 ( @{$omniscient->{'level2'}{$primary_tag_l2}{$id_l1}}) {
						my $level2_ID  = lc($feature_l2->_tag_value('ID'));
 						
 						if( $feature_l2->has_tag($original_attribute) and ! $feature_l2->has_tag('transcript_id') ){
							create_or_replace_tag($feature_l2, "transcript_id", $feature_l2->_tag_value($original_attribute) );
						}
					}
				}
			}
		}
	}
}

# @Purpose: add gene_id to all feature l1 if does not exist.
# @input: 3 =>  omniscient Hash reference, String = attribute to use to create gene_id
# @output none =>  The hash itself is modified
sub add_gene_id{
	my ($omniscient, $original_attribute) = @_;

	foreach my $primary_tag_l1 (keys %{$omniscient->{'level1'}}){ # primary_tag_l1 = gene or repeat etc...
		foreach my $id_l1 (keys %{$omniscient->{'level1'}{$primary_tag_l1}}){
			my $feature_l1  = $omniscient->{'level1'}{$primary_tag_l1}{$id_l1};
 						
 			if( $feature_l1->has_tag($original_attribute) and ! $feature_l1->has_tag('gene_id') ){
				create_or_replace_tag($feature_l1, "gene_id", $feature_l1->_tag_value($original_attribute) );
			}
		}
	}
}

#				   +------------------------------------------------------+
#				   |+----------------------------------------------------+|
#				   || 			FEATURES LOCATIONSATION					 ||
#				   |+----------------------------------------------------+|
#				   +------------------------------------------------------+

# looking the end and the start, the method check if two location overlap.
#A location os [Id, position1, position2]
# return t1 is location overlap
sub location_overlap{
	my($location1, $location2)=@_;
	my $overlap = undef;

	if (($location1->[1] <= $location2->[2]) and ($location1->[2] >= $location2->[1])){
		$overlap = 1;
	}

	return $overlap;
}

# looking the end and the start, the method check if two location overlap.
#A location os [Id, position1, position2]
# return the intersect of locations
sub location_overlap_update{
	my($location1, $location2)=@_;
	my $location = $location1;
	my $overlap = undef;

	if (($location1->[1] <= $location2->[2]) and ($location1->[2] >= $location2->[1])){
		$overlap = 1;
		if($location2->[1] < $location1->[1]){
			$location->[1] = $location2->[1]
		}
		if($location2->[2] > $location1->[2]){
			$location->[2] = $location2->[2]
		}
	}

	return $location, $overlap;
}

#Check if two genes have at least one L2 isoform which overlap at cds level.
# if no CDS we check i overlap at any other l3 feature.
sub check_gene_overlap_at_CDSthenEXON{
  my  ($hash_omniscient, $hash_omniscient2, $gene_id, $gene_id2)=@_;
  my $resu=undef;

	foreach my $l2_type (keys %{$hash_omniscient->{'level2'}} ){

		#check full CDS for each mRNA
		if(exists_keys($hash_omniscient,('level2', $l2_type, lc($gene_id)))){
			foreach my $mrna_feature (@{$hash_omniscient->{'level2'}{$l2_type}{lc($gene_id)}}){
				my $mrna_id1 = $mrna_feature->_tag_value('ID');

				if(exists_keys($hash_omniscient2,('level2', $l2_type, lc($gene_id2)))){
				    
				    foreach my $mrna_feature2 (@{$hash_omniscient2->{'level2'}{$l2_type}{lc($gene_id2)}}){ # from here bothe feature level2 are the same type
						my $mrna_id2 = $mrna_feature2->_tag_value('ID');
				   
					    #check all cds pieces
					    if(exists_keys($hash_omniscient,('level3', 'cds', lc($mrna_id1)))){
					      	if(exists_keys($hash_omniscient2,('level3', 'cds', lc($mrna_id2)))){
							    foreach my $cds_feature1 (@{$hash_omniscient->{'level3'}{'cds'}{lc($mrna_id1)}}){
							        foreach my $cds_feature2 (@{$hash_omniscient2->{'level3'}{'cds'}{lc($mrna_id2)}}){
							          
							        	if(($cds_feature2->start <= $cds_feature1->end) and ($cds_feature2->end >= $cds_feature1->start )){ # they overlap
							            	$resu="yes";last;
							          	}
							        }
							        if($resu){last;}
							    }

					      		if($resu){last;}
					      	}
					    }
					    elsif(! exists_keys($hash_omniscient2,('level3', 'cds', lc($mrna_id2)))){ # No CDS at all, check at exon / match level and if same level2 type
					    	
					    	foreach my $tag_l3 (keys %{$hash_omniscient->{'level3'}}){
					    		
					    		if(exists_keys($hash_omniscient,('level3', $tag_l3, lc($mrna_id1)))){					    			
					    			foreach my $feature1 (@{$hash_omniscient->{'level3'}{$tag_l3}{lc($mrna_id1)}}){
										
										if(exists_keys($hash_omniscient2,('level3', $tag_l3, lc($mrna_id2)))){					    					
					    					foreach my $feature2 (@{$hash_omniscient2->{'level3'}{$tag_l3}{lc($mrna_id2)}}){
					    						
					    						if(($feature2->start <= $feature1->end) and ($feature2->end >= $feature1->start )){ # they overlap
							            			$resu="yes";last;
							          			}
							          		}
								          	if($resu){last;}
								        }
								    }
								    if($resu){last;}
								}
							}

							if($resu){last;}
					    }
				    }

				    if($resu){last;}  
				}
			}

			if($resu){last;}
		}
	}
  return $resu;
}

# Check the start and end of level1 feature based on all features level2;
#return 1 if something modified
sub check_level1_positions {
	my ($hash_omniscient, $feature_l1, $verbose) = @_;
	my $result=undef;
	if(! $verbose){$verbose=0;}

	my $extrem_start=1000000000000;
	my $extrem_end=0;
	my $check_existence_feature_l2=undef;
	my $id_l1 = lc($feature_l1->_tag_value('ID'));

	foreach my $tag_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
    	if ( exists_keys ($hash_omniscient, ('level2', $tag_level2, $id_l1) ) ){
    		$check_existence_feature_l2=1;

	    	my $extrem_start_A=1000000000000;
			my $extrem_end_A=0;
	   		foreach my $feature ( @{$hash_omniscient->{'level2'}{$tag_level2}{$id_l1}}) {
	      		my $start=$feature->start();
	      		my $end=$feature->end();
	      		if ($start < $extrem_start_A){
	       			$extrem_start_A=$start;
	      		}
	      		if($end > $extrem_end_A){
	        		$extrem_end_A=$end;
	      		}
	      	}

	    	if ($extrem_start_A < $extrem_start){
	    		$extrem_start=$extrem_start_A;
	    	}
	    	if($extrem_end_A > $extrem_end){
	    		$extrem_end=$extrem_end_A;
	    	}
	    }
    }
    if(! $check_existence_feature_l2){
    	warn "WARNING check_level1_positions: NO level2 feature to check positions of the level1 feature ! @\n";
    }
    else{
	    # modify START if needed
	    if($feature_l1->start != $extrem_start){
	    	$feature_l1->start($extrem_start);
	    	$result=1;
	    	print "We modified the L1 LEFT extremity for the sanity the biological data!\n" if($verbose >= 3);
	    }

	    # modify END if needed
	    if($feature_l1->end != $extrem_end){
	    	$feature_l1->end($extrem_end);
	    	$result=1;
	    	print "We modified the L1 RIGHT extremity for the sanity the biological data!\n" if($verbose >= 3);
	    }
	}
	return $result;
}

# Check the start and end of level2 feature based on all features level3;
sub check_level2_positions {
	my ($hash_omniscient, $level2_feature)=@_;

	 my @values = $level2_feature->get_tag_values('ID');
     my $level2_feature_name = lc(shift @values) ;

	my $extrem_start=1000000000000;
	my $extrem_end=0;
	foreach my $tag_level3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
    	my $extrem_start_A=1000000000000;
		my $extrem_end_A=0;
		if( exists_keys ($hash_omniscient, ('level3', $tag_level3, $level2_feature_name) ) ){
	  		foreach my $feature ( @{$hash_omniscient->{'level3'}{$tag_level3}{$level2_feature_name}}) {
	   			my $start=$feature->start();
	   			my $end=$feature->end();
	   			if ($start < $extrem_start_A){
	   				$extrem_start_A=$start;
	   			}
	   			if($end > $extrem_end_A){
	      			$extrem_end_A=$end;
	   			}
	   		}
	   	}
   		if ($extrem_start_A < $extrem_start){
    		$extrem_start=$extrem_start_A;
    	}
   		if($extrem_end_A > $extrem_end){
    		$extrem_end=$extrem_end_A;
    	}
    }

    # modify START if needed
    if($level2_feature->start != $extrem_start){
    	$level2_feature->start($extrem_start);
    }

    # modify END if needed
    if($level2_feature->end != $extrem_end){
    	$level2_feature->end($extrem_end);
    }
}

#calcul the overlaping percentage betwwen 2 CDS list or 2 exon list etc...
# /!\ Be careful if you test the output, a overlaping gene can have a percentage overlap to 0. And if you test "if(featuresList_overlap)" and you have a 0, it will fail. So you have to check defined(featuresList_overlap)
sub featuresList_overlap {

	my ($listCDS1ref, $listCDS2ref)=@_;
	my $resu;

	###
	# sort the list
	my @listCDS1 = sort {$a->start <=> $b->start} @{$listCDS1ref};
	my @listCDS2 = sort {$a->start <=> $b->start} @{$listCDS2ref};
	# foreach my $t (@listCDS1){
	# 	print "list1: ".$t->start." ".$t->end."\n";
	# }
	# foreach my $t (@listCDS2){
	# 	print "list2: ".$t->start." ".$t->end."\n";
	# }
	my $size_overlap=0;
	my $cds1_size=0;
	foreach my $cds1 (@listCDS1){

		$cds1_size=$cds1_size+($cds1->end - $cds1->start)+1;
		my $starto;
		my $endo;

		foreach my $cds2 (@listCDS2){
			if($cds2->start > $cds1->end){ #we are after of the investigated cds.
				last;
			}
			elsif($cds2->end < $cds1->start){ # we are before investigated cds.
				next;
			}
			else{ #we are overlaping
				#check start
				if($cds1->start >= $cds2->start){
					$starto=$cds1->start;
				}
				else{$starto=$cds2->start;}
				#check end
				if($cds1->end >= $cds2->end){
					$endo=$cds2->end;
				}
				else{$endo=$cds1->end;}

				#calcul overlap;
				$size_overlap=$size_overlap+($endo - $starto + 1);
			}
		}
	}

	#Now calcul percentage overlap
	if($size_overlap != 0){
		$resu=($size_overlap*100)/$cds1_size;
		$resu = sprintf('%.0f', $resu);
		return $resu;
	}
	else{return undef;}
}


sub featuresList_identik {
	my ($list1, $list2)=@_;

	my @slist1 = sort {$a->start <=> $b->start} @{$list1};
	my @slist2 = sort {$a->start <=> $b->start} @{$list2};
	my $identik="true";

	if($#slist1 == $#slist2){
		my $cpt=0;

		while ($cpt <= $#slist1){

			my $feature1=$slist1[$cpt];
			my $feature2=$slist2[$cpt];
			#print $feature1->start." != ".$feature2->start." or  ".$feature1->end." != ".$feature2->end." or  ".$feature1->strand." ne ".$feature2->strand." or  ".$feature1->seq_id." ne ".$feature2->seq_id."\n";
			if( ($feature1->start != $feature2->start) or ($feature1->end != $feature2->end) or ($feature1->strand ne $feature2->strand) or ($feature1->seq_id ne $feature2->seq_id)){
				$identik=undef;last;
			}
			$cpt++;
		}
	}
	else{$identik=undef;}
	return $identik;
}

# @Purpose: Check the start and end of l1 l2 features from l3 features.
# @input: 2 => hash(omniscient hash), string(gene identifier)
# @output: none
sub check_record_positions {
	my ($hash_omniscient, $gene_id_raw)=@_;
	  
	my $gene_id = lc($gene_id_raw);
	my $ExtremStart=1000000000000;
	my $ExtremEnd=0;

  	foreach my $primary_tag_l2 (keys %{$hash_omniscient->{'level2'}} ){ 
	  	if (exists_keys($hash_omniscient, ('level2', $primary_tag_l2, $gene_id ) ) ){ 
		    foreach my $mrna_feature ( @{$hash_omniscient->{'level2'}{$primary_tag_l2}{$gene_id}} ) {
		      	my $l2_id = lc($mrna_feature->_tag_value('ID'));
		      	my $l2_ExtremStart=1000000000000;
	  			my $l2_ExtremEnd=0;
		      	foreach my $tag_l3 ( keys %{$hash_omniscient->{'level3'}} ) { 
	  				if ( exists_keys ( $hash_omniscient, ('level3', $tag_l3, $l2_id ) ) ){ 
		    			foreach my $feature_l3 ( @{$hash_omniscient->{'level3'}{$tag_l3}{$l2_id}} ) {

						    if ($feature_l3->start() < $l2_ExtremStart){
						       $l2_ExtremStart = $feature_l3->start();
						    }
						    if($feature_l3->end() > $l2_ExtremEnd){
						       $l2_ExtremEnd = $feature_l3->end();
			      			}
			      		}
		      		}
		      	}
		      	if ($mrna_feature->start != $l2_ExtremStart and $l2_ExtremStart != 1000000000000){
			      $mrna_feature->start($l2_ExtremStart);
			   	}
			  	if($mrna_feature->end != $l2_ExtremEnd and $l2_ExtremEnd != 0){
			   	 $mrna_feature->end($l2_ExtremEnd);
			  	}
			  	if ( $l2_ExtremStart < $ExtremStart ){
			  		$ExtremStart = $l2_ExtremStart;
			  	}
			  	if ($l2_ExtremEnd > $ExtremEnd  ){
			  		$ExtremEnd = $l2_ExtremEnd;
			  	}
		    }
		}


	  	my $gene_feature=$hash_omniscient->{'level1'}{'gene'}{$gene_id};
	  	if ($gene_feature->start != $ExtremStart and $ExtremStart != 1000000000000){
	      $gene_feature->start($ExtremStart);
	   	}
	  	if($gene_feature->end != $ExtremEnd and $ExtremEnd != 0){
	   	 $gene_feature->end($ExtremEnd);
	  	}
	}
}


# @Purpose: Check the start and end of gene feature based on its mRNAs and eventualy fix it.
# @input: 2 => hash(omniscient hash), string(gene identifier)
# @output: none
sub check_gene_positions {

  my ($hash_omniscient, $gene_id)=@_;

  #####
  #Modify gene start-end (have to check size of each mRNA)
  my $geneExtremStart=1000000000000;
  my $geneExtremEnd=0;
  foreach my $primary_tag_l2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
  	if (exists_keys($hash_omniscient, ('level2', $primary_tag_l2, lc($gene_id) ) ) ){ # check if they have mRNA avoiding autovivifcation
	    foreach my $mrna_feature ( @{$hash_omniscient->{'level2'}{$primary_tag_l2}{lc($gene_id)}}) {
	      	my $start=$mrna_feature->start();
	      	my $end=$mrna_feature->end();

	      if ($start < $geneExtremStart){
	        $geneExtremStart=$start;
	      }
	      if($end > $geneExtremEnd){
	        $geneExtremEnd=$end;
	      }
	    }
	}
  }
  my $gene_feature=$hash_omniscient->{'level1'}{'gene'}{lc($gene_id)};
  if ($gene_feature->start != $geneExtremStart){
      $gene_feature->start($geneExtremStart);
   }
  if($gene_feature->end != $geneExtremEnd){
    $gene_feature->end($geneExtremEnd);
  }
}

# From an omniscient, create a hash{primary_tag}{position} of feature. It is sort by seq id because position = contig1+
sub sort_by_seq_id {
	my ($omniscient)=@_;

	my %hash_sortBySeq;
	foreach my $tag_level1 (keys %{$omniscient->{'level1'}}){
	  foreach my $level1_id (keys %{$omniscient->{'level1'}{$tag_level1}}){

	    if (exists_keys($omniscient, ('level2','mrna',$level1_id)) ){ # check if they have mRNA avoiding autovivifcation
	      my $mrna_id  = lc($omniscient->{'level2'}{'mrna'}{$level1_id}[0]->_tag_value('ID'));

	      if (exists_keys($omniscient, ('level3', 'cds', $mrna_id)) ){ # check if they have cds avoiding autovivification. Allow to skip tRNA.
	        my $position=$omniscient->{'level1'}{$tag_level1}{$level1_id}->seq_id."".$omniscient->{'level1'}{$tag_level1}{$level1_id}->strand;
	        push (@{$hash_sortBySeq{$tag_level1}{$position}}, $omniscient->{'level1'}{$tag_level1}{$level1_id});
	      }
	    }
	  }
	}
	return \%hash_sortBySeq;
}

# @Purpose: Check if two genes have at least one mRNA isoform which overlap at cds level.
# @input: 4 => hash(omniscient), string(gene identifier), hash(omniscient), string(gene identifier)
# @output: 1 => undef || string(yes)
sub _two_features_overlap_two_hashes{
  my  ($hash1, $gene_id1, $hash2, $gene_id2)=@_;
  my $resu=undef;

  #check full CDS for each mRNA
  foreach my $mrna_feature (@{$hash1->{'level2'}{'mrna'}{lc($gene_id1)}}){
    foreach my $mrna_feature2 (@{$hash2->{'level2'}{'mrna'}{lc($gene_id2)}}){

      my $mrna_id1 = $mrna_feature->_tag_value('ID');
      my $mrna_id2 = $mrna_feature2->_tag_value('ID');

      #check all cds pieces
      foreach my $cds_feature1 (@{$hash1->{'level3'}{'cds'}{lc($mrna_id1)}}){
        foreach my $cds_feature2 (@{$hash2->{'level3'}{'cds'}{lc($mrna_id2)}}){

          if(($cds_feature2->start <= $cds_feature1->end) and ($cds_feature2->end >= $cds_feature1->start )){ # they overlap
            $resu="yes";last;
          }
        }
        if($resu){last;}
      }
      if($resu){last;}
    }
    if($resu){last;}
  }
  return $resu;
}

# @Purpose: The hash of reference will be the Hash target (HashT). The nkept name will come from the hash of reference.
# When an overlap is found, the ID/parent are fixed and we return 1 as a success !
# @input: 4 => object(gene feature), hash(omniscient), hash(omniscient), hash(sortBySeq)
# @output: 1 => undef || integer(1)
sub find_overlap_between_geneFeature_and_sortBySeqId {
	my ($geneFeature, $hash_source, $hashT, $hashT_sortBySeq )=@_;

	my $tag = $geneFeature->primary_tag;
	my $seqid = $geneFeature->seq_id;
	my $strand = $geneFeature->strand;
	my $gene_idS = $geneFeature->_tag_value('ID');

	#find overlap
	my $total_overlap=0;
	my $nb_feat_overlap=0;
	my @ListOverlapingGene=();
    foreach my $gene_featureT ( @{$hashT_sortBySeq->{$tag}{"$seqid$strand"}}){

    	my $gene_idT = $gene_featureT->_tag_value('ID');


    	if($gene_idT eq $gene_idS){ next;} # avoid to compare same feature if we are checking same omniscient

    	my ($start1,$end1) = get_longest_cds_start_end($hashT,$gene_idT); # look at CDS because we want only ioverlapinng CDS
    	my ($start2,$end2) = get_longest_cds_start_end($hash_source,$gene_idS); # look at CDS becaus ewe want only ioverlapinng CDS

    	if( ($start2 <= $end1) and ($end2 >= $start1) ){ #feature overlap considering extrem start and extrem stop. It's just to optimise the next step. Avoid to do the next step every time. So at the end, that test (current one) could be removed

            #now check at each CDS feature independently
            if (_two_features_overlap_two_hashes($hash_source,$gene_idS, $hashT, $gene_idT)){
              #print "These two features overlap without same id ! :\n".$geneFeature->gff_string."\n".$gene_featureT->gff_string."\n";
              $nb_feat_overlap++;

              push(@ListOverlapingGene, $gene_featureT);
            }
        }
    }

     # Now manage name if some feature overlap
    if( $nb_feat_overlap > 0){
    	my $reference_feature = shift(@ListOverlapingGene);
      	push(@ListOverlapingGene, $geneFeature);
        #print "$nb_feat_overlap overlapping feature found ! We will treat them now:\n";
        #print "We decided to keep that one: ".$reference_feature->gff_string."\n";

        my $gene_id_ref  = $reference_feature->_tag_value('ID');

        #change level2 parent for feature of level2 that have a feature of level1 in $ListToRemove list
        foreach my $featureToRemove (@ListOverlapingGene){

        	my $gene_id_to_remove  = lc($featureToRemove->_tag_value('ID'));

        	#######
        	#which hash the feature come from ?
        	my $currentHash=undef;
        	foreach my $tag_l1 (keys %{$hash_source->{'level1'}} ){ # primary_tag_key_level1 = gene or repeat etc...
    			if($hash_source->{'level1'}{$tag_l1}{$gene_id_to_remove} ){
    				$currentHash = $hash_source;
    			}
    		}
    		if(! $currentHash){$currentHash = $hashT;}
    		# ok now hash is choosen
    		################

	        foreach my $tag_level2 (keys %{$currentHash->{'level2'}}){

            	if (exists_keys($currentHash, ('level2',$tag_level2,$gene_id_to_remove)) ){ # check if they have cds avoiding autovivification.

	              	my @list_tmp_features = @{$currentHash->{'level2'}{$tag_level2}{$gene_id_to_remove}}; # As we will remove element of the list we cannot loop over it directly, we have to save the list in a temporary list;
			        foreach my $level2_feature (@list_tmp_features){ #replace Parent of each feature level2 by the new level1 reference
			            # Change parent feature
			            create_or_replace_tag($level2_feature,'Parent',$gene_id_ref);

			            #add it in other list
			            push (@{$currentHash->{'level2'}{$tag_level2}{lc($gene_id_ref)}},$level2_feature);

			            #remove mRNA from list <= not mandatory
			            my $mrna_id_to_remove = $level2_feature->_tag_value('ID');
		                my @tag_list=('all');
		                my @id_list=($gene_id_to_remove);my @id_list2=($mrna_id_to_remove);

			            remove_element_from_omniscient(\@id_list, \@id_list2, $currentHash, 'level2', 'false', \@tag_list);
              		}
            	}
          	}

          	foreach my $tag_level1 (keys %{$currentHash->{'level1'}}){ # remove the old feature level1 now
              	delete $currentHash->{'level1'}{$tag_level1}{$gene_id_to_remove}; # delete level1
          	}
        } #END FEATURE TO HANDLE
        ###
        # check end and start of the new feature
        my $gene_id=lc($reference_feature->_tag_value('ID'));
        check_gene_positions($hashT, $gene_id);
      	return 1;
    }
    else{return undef;}
}

# looking the end and the start, the method check if two features overlap.
sub check_if_feature_overlap{
	my($feature1, $feature2)=@_;
	my $result=undef;
	if (($feature1->start <= $feature2->end) and ($feature1->end >= $feature2->start)){
		$result="true";
	}

return $result
}

# Sort by locusID !!!!
# L1 => LocusID->level->typeFeature->ID =[ID,start,end]
# L2 and L3 => LocusID->level->typeFeature->Parent->ID = [ID,start,end]
#
#
sub sort_by_seq{
	my ($omniscient) = @_;

	my %hash_sortBySeq;
  
  	foreach my $tag_level1 (keys %{$omniscient->{'level1'}}){
    	foreach my $level1_id (keys %{$omniscient->{'level1'}{$tag_level1}}){
	    	my $level1_feature = $omniscient->{'level1'}{$tag_level1}{$level1_id};
	    	my $ID = $level1_feature->_tag_value('ID');
	    	my $strand="+";
	    	if($level1_feature->strand != 1){$strand = "-";}
	    	my $position_l1=$level1_feature->seq_id."".$strand;

	    	$hash_sortBySeq{$position_l1}{"level1"}{$tag_level1}{$level1_id} = [$ID, int($level1_feature->start), int($level1_feature->end)];
        }
	}
	return \%hash_sortBySeq;
}

#Sorting mixed strings => Sorting alphabetically first, then numerically
# how to use: my @y = sort by_number @x;
sub alphaNum {
    my ( $alet , $anum ) = $a =~ /([^\d]+)(\d+)/;
    my ( $blet , $bnum ) = $b =~ /([^\d]+)(\d+)/;
    ( $alet || "a" ) cmp ( $blet || "a" ) or ( $anum || 0 ) <=> ( $bnum || 0 )
}
1;
