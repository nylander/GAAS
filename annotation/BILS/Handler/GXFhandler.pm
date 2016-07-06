#!/usr/bin/perl -w

package BILS::Handler::GXFhandler ;
#use Time::HiRes;
use strict;
use Data::Dumper;
use Exporter qw(import);
use URI::Escape;
use Bio::Seq;

our $VERSION     = 1.00;
our @ISA         = qw(Exporter);
our @EXPORT_OK   = qw(nb_feature_level1 check_gene_positions find_overlap_between_geneFeature_and_sortBySeqId sort_by_seq_id webapollo_compliant extract_cds_sequence group_l1features_from_omniscient create_omniscient_from_idlevel2list get_feature_l2_from_id_l2_l1 print_duplicates remove_omniscient_elements_from_level2_feature_list featuresList_identik fill_omniscient_from_other_omniscient group_features_from_omniscient featuresList_overlap check_level1_positions check_level2_positions info_omniscient fil_cds_frame exists_keys gtf2gff_features_in_omniscient_from_level1_id_list _group_features_by_transcript_and_seq_id reconstruct_locus_without_transcripts_with_seq_id remove_element_from_omniscient append_omniscient merge_omniscients remove_omniscient_elements_from_level1_id_list fill_omniscient_from_other_omniscient_level1_id print_omniscient_from_level1_id_list check_if_feature_overlap remove_tuple_from_omniscient print_ref_list_feature print_omniscient create_or_replace_tag);
our %EXPORT_TAGS = ( DEFAULT => [qw()],
                 Ok    => [qw(nb_feature_level1 check_gene_positions find_overlap_between_geneFeature_and_sortBySeqId sort_by_seq_id webapollo_compliant extract_cds_sequence group_l1features_from_omniscient create_omniscient_from_idlevel2list get_feature_l2_from_id_l2_l1 print_duplicates remove_omniscient_elements_from_level2_feature_list featuresList_identik fill_omniscient_from_other_omniscient group_features_from_omniscient featuresList_overlap check_level1_positions check_level2_positions info_omniscient fil_cds_frame exists_keys gtf2gff_features_in_omniscient_from_level1_id_list _group_features_by_transcript_and_seq_id reconstruct_locus_without_transcripts_with_seq_id remove_element_from_omniscient append_omniscient merge_omniscients remove_omniscient_elements_from_level1_id_list fill_omniscient_from_other_omniscient_level1_id print_omniscient_from_level1_id_list check_if_feature_overlap remove_tuple_from_omniscient print_ref_list_feature print_omniscient create_or_replace_tag)]);
=head1 SYNOPSIS



=head1 DESCRIPTION

	A library to convert
	Inherits from

	Dont take in account repeat and multi parent feature!!!

=cut


###################
#
# Print methods
#
###################
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

	#################
	# == LEVEL 1 == #
	#################
	foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
		foreach my $id_tag_key_level1 (keys %{$hash_omniscient->{'level1'}{$primary_tag_key_level1}}){
			$gffout->write_feature($hash_omniscient->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1}); # print feature

			#################
			# == LEVEL 2 == #
			#################
			foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

				if ( exists ($hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1} ) ){
					foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}) {
						$gffout->write_feature($feature_level2);

						#################
						# == LEVEL 3 == #
						#################
						my @temp = $feature_level2->get_tag_values('ID');
						my $level2_ID = lc(shift @temp);

						######
						# FIRST EXON
						if ( exists ($hash_omniscient->{'level3'}{'exon'}{$level2_ID} ) ){
							foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{'exon'}{$level2_ID}}) {
								$gffout->write_feature($feature_level3);
							}
						}
						###########
						# SECOND CDS
						if ( exists ($hash_omniscient->{'level3'}{'cds'}{$level2_ID} ) ){
							foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{'cds'}{$level2_ID}}) {
								$gffout->write_feature($feature_level3);
							}
						}

						############
						# THEN ALL THE REST
						foreach my $primary_tag_key_level3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
							if (($primary_tag_key_level3 ne 'cds') and ($primary_tag_key_level3 ne 'exon')) {
								if ( exists ($hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID} ) ){
									foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID}}) {
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

##################
#
# END Print methods
#
###################

# omniscient is a hash containing a whole gXf file in memory sorted in a specific way (3 levels)
sub gtf2gff_features_in_omniscient_from_level1_id_list {

	my ($hash_omniscient, $level_id_list, $gffout) = @_  ;

	#################
	# == LEVEL 1 == #
	#################
	foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
		foreach my $id_tag_key_level1 (@$level_id_list){

			if(exists ($hash_omniscient->{'level1'}{$primary_tag_key_level1}{lc($id_tag_key_level1)})){
				my $feature_level1 = $hash_omniscient->{'level1'}{$primary_tag_key_level1}{lc($id_tag_key_level1)};
				if(! $feature_level1->has_tag('ID')){

					if($feature_level1->has_tag('gene_id')){
						my @temp = $feature_level1->get_tag_values('gene_id');
						my $level1_ID = shift @temp;
						$feature_level1->add_tag_value('ID',$level1_ID)
					}
					else{
						warn "Cannot retrieve the id feature of the following feature: ".gff_string($feature_level1);
					}
				}

				#################
				# == LEVEL 2 == #
				#################
				foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

					if ( exists ($hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1} ) ){
						foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}) {

							##manage primary_tag
							if(lc($feature_level2->primary_tag) eq "transcript"){
								$feature_level2->primary_tag('mrna');
							}

							##manage ID
							my $level2_ID;
							if($feature_level2->has_tag('ID')){
								my @temp = $feature_level2->get_tag_values('ID');
								$level2_ID = shift @temp;
							}
							if(! $feature_level2->has_tag('ID')){
								if($feature_level2->has_tag('transcript_id')){
									my @temp = $feature_level2->get_tag_values('transcript_id');
									$level2_ID = shift @temp;
									$feature_level2->add_tag_value('ID',$level2_ID)
								}
								else{
									warn "Cannot retrieve the id feature of the following feature: ".gff_string($feature_level2);
								}
							}
							## Manage Parent
							if(! $feature_level2->has_tag('Parent')){

								if($feature_level2->has_tag('gene_id')){
									my @temp = $feature_level2->get_tag_values('gene_id');
									my $level2_Parent = shift @temp;
									$feature_level2->add_tag_value('Parent',$level2_Parent)
								}
								else{
									warn "Cannot retrieve the parent feature of the following feature: ".gff_string($feature_level2);
								}
							}

							#################
							# == LEVEL 3 == #
							#################
							foreach my $primary_tag_key_level3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
								my $feature_counter=1;
								if ( exists ($hash_omniscient->{'level3'}{$primary_tag_key_level3}{lc($level2_ID)} ) ){
									foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{$primary_tag_key_level3}{lc($level2_ID)}}) {

										# Manage Parent Name
										if(! $feature_level3->has_tag('Parent')){

											if($feature_level3->has_tag('transcript_id')){
												my @temp = $feature_level3->get_tag_values('transcript_id');
												my $level3_Parent = shift @temp;
												$feature_level3->add_tag_value('Parent',$level3_Parent)
											}
											else{
												warn "Cannot retrieve the parent feature of the following feature: ".gff_string($feature_level2);
											}
										}
										# Manage ID
										if(! $feature_level3->has_tag('ID')){
											if (lc($primary_tag_key_level3) eq "exon"){

												if($feature_level3->has_tag('exon_id')){
													my @temp = $feature_level3->get_tag_values('exon_id');
													my $level3_ID = shift @temp;
													$feature_level3->add_tag_value('ID',$level3_ID)
												}
												else{

													if($feature_level3->has_tag('transcript_id')){
														my @temp = $feature_level3->get_tag_values('transcript_id');
														my $level3_ID = shift @temp;
														$feature_level3->add_tag_value('ID',$level3_ID."-". $feature_level3->primary_tag . "-" . $feature_counter );
														$feature_counter++;
													}
													else{
														warn "Cannot retrieve the ID feature of the following feature: ".gff_string($feature_level3);
													}
												}
											}
											elsif(lc($primary_tag_key_level3) eq "cds"){

												if($feature_level3->has_tag('protein_id')){
													my @temp = $feature_level3->get_tag_values('protein_id');
													my $level3_ID = shift @temp;
													$feature_level3->add_tag_value('ID',$level3_ID)
												}
												else{

													if($feature_level3->has_tag('transcript_id')){
														my @temp = $feature_level3->get_tag_values('transcript_id');
														my $level3_ID = shift @temp;
														$feature_level3->add_tag_value('ID',$level3_ID."-". $feature_level3->primary_tag . "-" . $feature_counter );
														$feature_counter++;
													}
													else{
														warn "Cannot retrieve the ID feature of the following feature: ".gff_string($feature_level3);
													}
												}
											}
											else{

												if($feature_level3->has_tag('transcript_id')){
													my @temp = $feature_level3->get_tag_values('transcript_id');
													my $level3_ID = shift @temp;
													$feature_level3->add_tag_value('ID',$level3_ID."-". $feature_level3->primary_tag . "-" . $feature_counter );
													$feature_counter++;
												}
												else{
													warn "Cannot retrieve the ID feature of the following feature: ".gff_string($feature_level3);
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
	}
}

# omniscient is a hash containing a whole gXf file in memory sorted in a specific way (3 levels)
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

# omniscient is a hash containing a whole gXf file in memory sorted in a specific way (3 levels)
# Merge two omniscients. Features are added if ID was not already existsing. Starting from level1 to level3. If level1 already exists we do not go to the next level.
sub fill_omniscient_from_other_omniscient {

	my ($hash_omniscient, $omniscient_to_append)=@_;

	#################
	# == LEVEL 1 == #
	#################
	foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
		foreach my $id_tag_key_level1 (keys %{$hash_omniscient->{'level1'}{$primary_tag_key_level1}}){
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
							my $level2_ID = lc( $feature_level2->_tag_value('ID') );

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

# omniscient is a hash containing a whole gXf file in memory sorted in a specific way (3 levels)
# put data from hash_omniscient1 in hash_omniscient2
# Features are added if ID was not already existsing. Starting from level1 to level3. Each level is filled independently.
sub merge_omniscients {

	my ($hash_omniscient1, $hash_omniscient2)=@_;

	#################
	# == LEVEL 1 == #
	#################
	foreach my $primary_tag_key_level1 (keys %{$hash_omniscient1->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
		foreach my $id_tag_key_level1 (keys %{$hash_omniscient1->{'level1'}{$primary_tag_key_level1}}){
			if ( ! exists ($hash_omniscient2->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1} ) ){
					$hash_omniscient2->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1} = $hash_omniscient1->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1}; # print feature
			}
			else{print "INFO level1:  $id_tag_key_level1 already exist in file 1\n";}
		}
	}

	#################
	# == LEVEL 2 == #
	#################
	foreach my $primary_tag_key_level2 (keys %{$hash_omniscient1->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
		foreach my $id_tag_key_level2 (keys %{$hash_omniscient1->{'level2'}{$primary_tag_key_level2}}){

			if (! exists ($hash_omniscient2->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level2} ) ){ #Non present in hash2, we create a list with one element
					@{$hash_omniscient2->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level2}} = @{$hash_omniscient1->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level2}};
			}
			else{ # header of the list exist, does all the feature exists ?

				print "INFO level2:  Parent $id_tag_key_level2 already exist in file 1. Should we add a new mRNA ?\n";
				foreach my $feature_level2_hash1 ( @{$hash_omniscient1->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level2}}) {
					my $level2_exists="no";

					my @temp = $feature_level2_hash1->get_tag_values('ID');
					my $level2_ID_hash1 = lc(shift @temp);

					# Check if feature exists in hash2
					foreach my $feature_level2_hash2 ( @{$hash_omniscient2->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level2}}) {
						my @temp = $feature_level2_hash2->get_tag_values('ID');
						my $level2_ID_hash2 = lc(shift @temp);
						if($level2_ID_hash2 eq $level2_ID_hash1){
							my $level2_exists="yes";
						}
					}
					# feature doesnt exist in list of hash2, so we add it
					if ($level2_exists eq "no"){
						push (@{$hash_omniscient2->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level2}}, $feature_level2_hash1);
					}
				}
			}
		}
	}

	#################
	# == LEVEL 3 == #
	#################
	foreach my $primary_tag_key_level3 (keys %{$hash_omniscient1->{'level3'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
		foreach my $id_tag_key_level3 (keys %{$hash_omniscient1->{'level3'}{$primary_tag_key_level3}}){
			if (! exists ($hash_omniscient2->{'level3'}{$primary_tag_key_level3}{$id_tag_key_level3} ) ){ #Non present in hash2, we create a list with one element
					@{$hash_omniscient2->{'level3'}{$primary_tag_key_level3}{$id_tag_key_level3}} = @{$hash_omniscient1->{'level3'}{$primary_tag_key_level3}{$id_tag_key_level3}};
			}
			else{
				print "INFO level3:  $id_tag_key_level3 already exist in file 1. Should we add a new level3 features ?\n";
				foreach my $feature_level3_hash1 ( @{$hash_omniscient1->{'level3'}{$primary_tag_key_level3}{$id_tag_key_level3}}) {
					my $level3_exists="no";

					my @temp = $feature_level3_hash1->get_tag_values('ID');
					my $level3_ID_hash1 = lc(shift @temp);

					# Check if feature exists in hash2
					foreach my $feature_level3_hash2 ( @{$hash_omniscient2->{'level3'}{$primary_tag_key_level3}{$id_tag_key_level3}}) {
						my @temp = $feature_level3_hash2->get_tag_values('ID');
						my $level3_ID_hash2 = lc(shift @temp);
						if($level3_ID_hash2 eq $level3_ID_hash1){
							my $level3_exists="yes";
						}
					}
					# feature doesnt exist in list of hash2, so we add it
					if ($level3_exists eq "no"){
						push (@{$hash_omniscient2->{'level3'}{$primary_tag_key_level3}{$id_tag_key_level3}}, $feature_level3_hash1);
					}
				}
			}
		}
	}
}

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
								my @temp = $feature_level2->get_tag_values('ID');
								my $level2_ID = lc(shift @temp);



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

# omniscient is a hash containing a whole gXf file in memory sorted in a specific way (3 levels)
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
						my $feature_Parent = lc($feature->_tag_value('Parent'));

						if($level2_ID eq $feature_ID){

							#################
							# == LEVEL 3 == #
							#################
							foreach my $primary_tag_key_level3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...
								if( exists_keys($hash_omniscient, ('level3', $primary_tag_key_level3, $level2_ID)) ){
									delete $hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID} # delete level3
								}
							}
						my @id_concern_list=($feature_Parent);
						my @id_list_to_remove=($feature_ID);
						my @list_tag_key=('all');
						remove_element_from_omniscient(\@id_concern_list, \@id_list_to_remove, $hash_omniscient, 'level2','false', \@list_tag_key);
						#delete $hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1} # delete level2
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
					if($id_to_remove eq $id_key){
						delete $hash_omniscient->{$level}{$tag_key}{$id_to_remove}; #REMOVE THAT KEY-VALUE pair
					}
				}
			}
		}
	}
}

# from omniscient: remove feature from "feature list" of level2 or level3 with id present in $id_list_to_remove
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

				foreach my $feature (@{$hash_omniscient->{$level}{$tag_key}{$id_concern}}){
					my @values = $feature->get_tag_values('ID');
					my $id = lc(shift @values);
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
					@{$hash_omniscient->{$level}{$tag_key}{$id_concern}}=@listok;
				}
			}
		}
	}
}

sub create_omniscient {

	my ($level1,$level2,$level3)=@_;

	my $omniscient;

	foreach my $feature (@$level1){
		my $id = lc($feature->_tag_value('ID'));
		$omniscient->{"level1"}{lc($feature->primary_tag)}{$id}=$feature;
	}
	foreach my $feature (@$level2){
		my $id = lc($feature->get_tag_values('Parent'));
		push(@{$omniscient->{"level2"}{lc($feature->primary_tag)}{$id}}, $feature);###
	}
	foreach my $feature (@$level3){
		my $id = lc( $feature->get_tag_values('Parent'));
		push(@{$omniscient->{"level3"}{lc($feature->primary_tag)}{$id}}, $feature);
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

sub append_omniscient {

	my ($omniscient, $level1,$level2,$level3)=@_;

	foreach my $feature (@$level1){
		my $primaryTag = lc($feature->primary_tag);
		my @values = $feature->get_tag_values('ID');
		my $id = lc(shift @values) ;
		if( ! exists_keys($omniscient, ('level1', $primaryTag, $id)) ){
			$omniscient->{"level1"}{$primaryTag}{$id}=$feature;
		}
	}
	foreach my $feature (@$level2){ # if exist, try to append the list
		my $primaryTag = lc($feature->primary_tag);
		my @values = $feature->get_tag_values('Parent');
		my $parent_id = lc(shift @values) ;
		if( ! exists_keys($omniscient, ('level2', $primaryTag, $parent_id)) ){
			push(@{$omniscient->{"level2"}{$primaryTag}{$parent_id}}, $feature);###
		}
		else{ # append element in the list if not existing
			my $exist_in_list="no";
			my @values = $feature->get_tag_values('ID');
			my $id = lc(shift @values) ;
			foreach my $feature_original (@{$omniscient->{"level2"}{$primaryTag}{$parent_id}}){
				my @original_values = $feature_original->get_tag_values('ID');
				my $original_id = lc(shift @values);
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
		my @values = $feature->get_tag_values('Parent');
		my $parent_id = lc(shift @values) ;
		if( ! exists_keys($omniscient, ('level3', $primaryTag, $parent_id)) ){
			push(@{$omniscient->{"level3"}{$primaryTag}{$parent_id}}, $feature);
		}
		else{ # append element in the list if not existing
			my $exist_in_list="no";
			my @values = $feature->get_tag_values('ID');
			my $id = lc(shift @values) ;
			foreach my $feature_original (@{$omniscient->{"level3"}{$primaryTag}{$parent_id}}){
				my @original_values = $feature_original->get_tag_values('ID');
				my $original_id = lc(shift @values);
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

#Transform omniscient data to be Webapollo compliant
sub webapollo_compliant {
		my ($hash_omniscient) = @_  ;

	#################
	# == LEVEL 1 == #
	#################
	foreach my $primary_tag_key_level1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
		foreach my $id_tag_key_level1 (keys %{$hash_omniscient->{'level1'}{$primary_tag_key_level1}}){
			webapollo_rendering($hash_omniscient->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1});

			#################
			# == LEVEL 2 == #
			#################
			foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

				if ( exists ($hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1} ) ){
					foreach my $feature_level2 ( @{$hash_omniscient->{'level2'}{$primary_tag_key_level2}{$id_tag_key_level1}}) {
						webapollo_rendering($feature_level2);

						#################
						# == LEVEL 3 == #
						#################
						my $level2_ID  = lc($feature_level2->_tag_value('ID'));

						foreach my $primary_tag_key_level3 (keys %{$hash_omniscient->{'level3'}}){ # primary_tag_key_level3 = cds or exon or start_codon or utr etc...

							if ( exists ($hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID} ) ){
								foreach my $feature_level3 ( @{$hash_omniscient->{'level3'}{$primary_tag_key_level3}{$level2_ID}}) {
									webapollo_rendering($feature_level3);
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
sub webapollo_rendering {

	my ($feature)=@_;

	## check primary tag
	my $primary_tag = lc($feature->primary_tag);

	if($primary_tag eq 'cds'){
		$feature->primary_tag('CDS');
	}
	if($primary_tag eq 'exon'){
		$feature->primary_tag('exon');
	}
	if($primary_tag eq 'three_prime_utr'){
		$feature->primary_tag('three_prime_UTR');
	}
	if($primary_tag eq 'five_prime_utr'){
		$feature->primary_tag('five_prime_UTR');
	}
	if($primary_tag eq 'utr'){
		$feature->primary_tag('UTR');
	}
	if($primary_tag eq 'mrna'){
		$feature->primary_tag('mRNA');
	}
	if($primary_tag eq 'gene'){
		$feature->primary_tag('gene');
	}

	## check attribute
	if($feature->has_tag('product')){
		my @values = $feature->get_tag_values('product');
		$feature->add_tag_value('description', @values);
		$feature->remove_tag('product');
	}
}

#Transform omniscient data to be Webapollo compliant
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

# looking the end and the start, the method check if two features overlap.
sub check_if_feature_overlap{
	my($feature1, $feature2)=@_;
	my $result=undef;
	if (($feature1->start <= $feature2->end) or ($feature1->end >= $feature2->start)){
		$result="true";
	}

return $result
}

sub create_or_replace_tag{

	my ($feature, $tag, $value)=@_;

	if ($feature->has_tag($tag) ) {
			$feature->remove_tag($tag);
        	$feature->add_tag_value($tag,$value);
	}
	else{
		$feature->add_tag_value($tag,$value);
	}
}

############################## Based on Marc idea =====

sub reconstruct_locus_without_transcripts_with_seq_id {
my 	@gene;
my 	@mRNA;
	# All features belonging to the same gene
	my ($features, $NewGeneID) = @_; #Features are

	# Take all features and bin them by transcript_id

	my %bin = _group_features_by_transcript_and_seq_id(@$features);

	# Build transcript from all features in the group

	my @answer = ();

	# Some annotations (like lift-overs) can be split across scaffolds/contigs
	my $is_split = 0;
	my $appendix = ""; # A variable to store appendices for split genes
	# Check wether this gene is split across scaffolds
	if ( (keys %bin) > 1) {
		$is_split = 1;
	}

	my $seq_counter = 0;
	while (my ($seq_id,$transcript_hash) = each %bin) {

		# !!! These are all transcripts on the same sequence !!!

		# Collect transcripts for gene reconstruction
		my @transcripts = ();
		my $exon_counter = 0;


		while (my ($transcript_id,$features) = each %$transcript_hash) {

			my ($t_feature, $NewGeneIDtmp) = _build_transcript_from_exons_with_seq_id($features, $NewGeneID);
			$NewGeneID=$NewGeneIDtmp; # keep track If Gene Id modified

			if ($is_split == 1) {
				my $seqid=$t_feature->seq_id();
				$appendix = "partial_part-$seqid" ;
				$t_feature = _update_feature_id($t_feature,$appendix) ;
			}

			push(@transcripts,$t_feature);
			unshift (@{ $bin{$seq_id}{$transcript_id} }, $t_feature ) ;
		}

		# Build the gene container
		my $gene_feature = _build_gene_from_transcripts_with_seq_id(@transcripts);

		# Gene into array
		push(@gene,$gene_feature);
		push(@mRNA,@transcripts);

	} # End sequence region

	return \@gene,\@mRNA,$NewGeneID;
}

sub _update_feature_id {

	my $feature = shift ;
	my $appendix = shift;

	my @tags = ( 'ID' , 'Parent' ) ;

	foreach my $tag (@tags) {

		if ($feature->has_tag($tag) ) {
			my @temp = $feature->get_tag_values($tag);
			my $current_id = shift @temp;
	        	my $new_id = $current_id . "_" . $appendix ;
			$feature->remove_tag($tag);
        		$feature->add_tag_value($tag,$new_id);
		}
	}

	return $feature;
}

# Deals with features across contigs
sub _group_features_by_transcript_and_seq_id {

	my @features = @_;

	my %bin ;

	foreach my $feature (@features) {

		if ($feature->start < 1) {
			$feature->start(1);
		}

		my $source_tag = $feature->source_tag;
		my $primary_tag = $feature->primary_tag;

		# If a gene is split across multiple sequences,
		# we need to build separate genes for each sequence -
		# this keeps track of that.

		my $seq_id = $feature->seq_id;

		# We need to reorganize feature relationships into parent/child:
		my $this_id = undef;
		my $parent_id = undef;
		my $transcript_id=undef;

		# We always need to know about the transcript_id, record it...
		if($feature->has_tag('Parent')){ #gff case
			$transcript_id = $feature->_tag_value('Parent');
		}
		elsif($feature->has_tag('transcript_id')){ #gtf case
			$transcript_id = $feature->_tag_value('transcript_id');
		}
		else{
			warn "No id found for this feature ".gff_string($feature);
		}


		# Group features by transcripts/mRNA and seq_id

		if ( !exists( $bin{$seq_id} ) ) {
			$bin{$seq_id} = {};
		}

		if ( !exists( $bin{$seq_id}{$transcript_id} ) ) {
			$bin{$seq_id}{$transcript_id} = [];
		}

		# Some annotation files (eg. lift-overs) may
		# included negative coordinates, need to fix
		push( @{ $bin{$seq_id}{$transcript_id} }, $feature );
	}

	return %bin ;

}

sub _build_transcript_from_exons_with_seq_id {

	my ($features, $NewGeneID) = @_;

	my $transcript_start = 0;
	my $transcript_end = 0;
	my $transcript_strand = undef;
	my $seq_id = undef;
	my $source_tag = undef;
	my $transcript_id = undef;
	my $primary_tag = undef;
	my $gene_id = undef;

	# All features belonging to the same transcript_id
	# Used to calculate transcript coordinates.
	foreach my $f (@{$features}) {

		if($f->has_tag('gene_id')){ #gtf case
			$gene_id = $f->_tag_value('gene_id');
		}
		elsif($f->has_tag('Parent')){ #gff case
			$gene_id = $f->_tag_value('Parent');
			print "GXF check what we should do in that case\n";
		}

		if($f->has_tag('transcript_id')){ #gtf case
			$transcript_id = $f->_tag_value('transcript_id');
		}
		elsif($f->has_tag('ID')){ #gff case
			$transcript_id = $f->_tag_value('ID');
			print "GXF check what we should do in that case2\n";

		}
		else{
			print "GXFhandler Warning: No transcript_id neither ID attribute. So we will create one.\n";
			print "actually it's not yet implemeted.... Sorry\n";
		}
		# CDS never contain more information than exons, skip
		#if ($f->primary_tag eq 'exon') {
			$transcript_strand = $f->strand;
			$seq_id = $f->seq_id;
			$source_tag = $f->source_tag;

			if ($transcript_start == 0) {
				$transcript_start = $f->start;
				$transcript_end = $f->end;
				$transcript_strand = $f->strand;
				$seq_id = $f->seq_id;
				$source_tag = $f->source_tag;
			}

			if ($transcript_start > $f->start) {
				$transcript_start = $f->start ;
			}

			if ($transcript_end < $f->end) {
				$transcript_end = $f->end ;
			}
		#}
	} # end transcript features


	## THIS NEEDS CHECKING!
	# We assume that EnsEMBL source_tags comply with SO terms - except 'protein_coding'

	#No gene_id ! We need one to create the mRNA
	if(!$gene_id){ # No information we must create a new gene_id
		print "GXFhandler Warning: No gene_id neither Parent attribute. So we will create one.\n";
		$NewGeneID++;
		$gene_id="gene".$NewGeneID;
	}

	if ($source_tag =~ /.*RNA.*/) {
		$primary_tag = $source_tag ;
	} else {
		$primary_tag = "mRNA" ;
	}

	my $t_feature = Bio::SeqFeature::Generic->new(-start => $transcript_start, -primary_tag => $primary_tag , -frame => '.' , -end => $transcript_end, -strand => $transcript_strand , -seq_id => $seq_id, -source_tag => $source_tag, -tag => { 'ID' => $transcript_id , 'Parent' => $gene_id }) ;

	return $t_feature, $NewGeneID;

}

sub _build_gene_from_transcripts_with_seq_id {

	my @transcripts = @_ ;

	# Define the gene container parameters:
	my $gene_start = 0;
	my $gene_end = 0;
	my $gene_strand = undef;
	my $gene_id = undef;
	my $contig = undef;
	my $source_tag = undef;
	my $seq_id = undef;

	foreach my $t_feature (@transcripts) {

		$source_tag = $t_feature->source_tag;
		$seq_id = $t_feature->seq_id;

		my @gvalues = $t_feature->get_tag_values('Parent');
		$gene_id = shift @gvalues ;

		# Re-size the gene container...
		if ($gene_start == 0) {
			$gene_start = $t_feature->start;
			$gene_end = $t_feature->end;
			$gene_strand = $t_feature->strand;
			$contig = $t_feature->seq_id;
		}
		if ($gene_start > $t_feature->start) {
			$gene_start = $t_feature->start;
		}

		if ($gene_end < $t_feature->end) {
			$gene_end = $t_feature->end;
		} # end gene resize

	}

	my $gene_feature = Bio::SeqFeature::Generic->new(-start => $gene_start, -primary_tag => 'gene' , -frame => '.', -end => $gene_end, -strand => $gene_strand , -seq_id => $contig, -source_tag => $source_tag, -tag => { "ID" => $gene_id }) ;

	return $gene_feature;

}

#============= END Marc idea develloped ===============

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

			        $resu{$tag}=$resu{$tag}+$nb;
     			}
    		}
  		}
	}
	foreach my $tag (keys %resu){
		print "There is $resu{$tag} $tag\n";
	}
}

#============
#check if reference exists in hash. Deep infinite : hash{a} or hash{a}{b} or hash{a}{b}{c}, etc.
# usage example: exists_keys($hash_omniscient,('level3','cds',$level2_ID)
sub _exists {
    my ($hash, @keys, @exists) = @_;

    foreach my $el (@keys){ #test all keys... no more !
    	if (exists $hash->{$keys[0]}) {

    	    my $key = shift @keys;
    	    push @exists, 1, _exists($hash->{$key}, @keys);
   		}
	}
    return @exists;
}

sub exists_keys {
    return &_exists == $#_ ? 1 : ();
}
#============

# Check the start and end of level1 feature based on all features level2;
sub check_level1_positions {
	my ($hash_omniscient, $level1_feature)=@_;

	my $level1_feature_name = lc( $level1_feature->_tag_value('ID'));

	my $extrem_start=1000000000000;
	my $extrem_end=0;
	my $check_existence_feature_l2=undef;

	foreach my $tag_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
    	if ( exists ($hash_omniscient->{'level2'}{$tag_level2}{$level1_feature_name} ) ){
    		$check_existence_feature_l2=1;

	    	my $extrem_start_A=1000000000000;
			my $extrem_end_A=0;
	   		foreach my $feature ( @{$hash_omniscient->{'level2'}{$tag_level2}{$level1_feature_name}}) {
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
    	warn "WARNING: NO level2 feature to check positions of the level1 feature ! @\n";
    }
    else{
	    # modify START if needed
	    if($level1_feature->start != $extrem_start){
	    	$level1_feature->start($extrem_start);
	    }

	    # modify END if needed
	    if($level1_feature->end != $extrem_end){
	    	$level1_feature->end($extrem_end);
	    }
	}
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
						my @temp = $feature_level2->get_tag_values('ID');
						my $level2_ID = lc(shift @temp);


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
			push(@{$group{$seq_id}}, $hash_omniscient->{'level1'}{$primary_tag_key_level1}{$id_tag_key_level1});

		}
	}
	return \%group;
}

# print duplicate hash
sub print_duplicates {
	my ($duplicate_omniscient, $hash_omniscient, $gffout) = @_  ;

	foreach my $level (keys %{$duplicate_omniscient}){ # primary_tag_key_level1 = gene or repeat etc...
		foreach my $primary_tag (keys %{$duplicate_omniscient->{$level}}){
			my $nb_by_pt=0;
			my $nb_feat_pt=0;
			foreach my $id (keys %{$duplicate_omniscient->{$level}{$primary_tag}}){
				$nb_feat_pt++;
				foreach my $feature (@{$duplicate_omniscient->{$level}{$primary_tag}{$id}}){
					$nb_by_pt++;
					$gffout->write_feature($feature); # print feature
				}
			}
			print "We found $nb_feat_pt duplicated $primary_tag feature for a total of $nb_by_pt duplicates.\n";
		}
	}
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

# @Purpose: from a omniscient and a gene_id, will get back the extrem value for start and end
# @input: 2 => hash(omniscient), string(gene identifier)
# @output: 2 => integer(extrem start position), integer(extrem end position)
sub get_longest_cds_start_end {
  my  ($hash_omniscient,$gene_id)=@_;
  my $resu_start=100000000000;
  my $resu_end=0;

  #check full CDS for each mRNA
  foreach my $mrna_feature (@{$hash_omniscient->{'level2'}{'mrna'}{lc($gene_id)}}){
    my @values = $mrna_feature->get_tag_values('ID');
    my $mrna_id = shift @values;
    my $extrem_start=100000000000;
    my $extrem_end=0;

    #check all cds pieces
    foreach my $cds_feature (@{$hash_omniscient->{'level3'}{'cds'}{lc($mrna_id)}}){
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
1;
