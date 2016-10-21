#!/usr/bin/perl -w

package BILS::Handler::GFF3handler ;

use strict;
use warnings;
use Data::Dumper;
use Clone 'clone';
use BILS::Handler::GXFhandler qw(:Ok);
use Exporter qw(import);

our $VERSION     = 1.00;
our @ISA         = qw(Exporter);
our @EXPORT_OK   = qw(check_mrna_positions modelate_utr_and_cds_features_from_exon_features_and_cds_start_stop create_omniscient slurp_gff3_file_JD create_omniscient_from_feature_list);
our %EXPORT_TAGS = ( DEFAULT => [qw()],
                 Ok    => [qw(check_mrna_positions modelate_utr_and_cds_features_from_exon_features_and_cds_start_stop create_omniscient slurp_gff3_file_JD create_omniscient_from_feature_list)]);
=head1 SYNOPSIS



=head1 DESCRIPTION

	A library to convert 
	Inherits from 

=cut	

# Save in omniscient hash (sorted in a specific way (3 levels)) a whole gff3 file
sub slurp_gff3_file_JD {
	
	my ($self, $file) = @_  ;
	
	my $gffio = Bio::Tools::GFF->new(-file => $file, -gff_version => 3);	

	my $start_run = time();

	### Handle to not print to much warning
	my %WARNS;
	my $nbWarnLimit=5;
  	local $SIG{__WARN__} = sub {
    my $message = shift;
    my @thematic=split /@/,$message ;
    
    $WARNS{$thematic[0]}++;
    	if ($WARNS{$thematic[0]} <= $nbWarnLimit){
			print $message;
		}
		if($WARNS{$thematic[0]} == $nbWarnLimit){
			print "$thematic[0] ************** Too much WARNING message we skip the next **************\n";
		}
  	}; 
  	########### END WARNING LOCAL METHOD


	my %mRNAGeneLink; #Hast that keep track about link between l2 and l1
	my %omniscient; #Hast where all the feature will be saved
	my %duplicate;# Hash to store any counter. Will be use to create a new uniq ID 
	my %miscCount;# Hash to store any counter. Will be use to create a new uniq ID

	#read every lines
	while( my $feature = $gffio->next_feature()) {
		manage_one_feature($feature, \%omniscient, \%mRNAGeneLink, \%duplicate, \%miscCount);
    }
    #close the file
    $gffio->close();

    #Check if duplicate detected:
    my $keyExist = keys %duplicate;
    if($keyExist){#print result
    	my $stringPrint .= "################################################################################\n# Achthung /\\ Attention /\\ Be carefull => Identical feature found !            #\n".
      	"# (Same chr/contig/scaffold, same position, same ID). They have been removed ! #\n################################################################################\n\n";
      	print $stringPrint;
      	my $gffout= Bio::Tools::GFF->new( -fh => \*STDOUT );
      	print_duplicates(\%duplicate, \%omniscient, $gffout);
      	print "We removed them from the regular output.\n";
    }

    #Check relationship between mRNA and gene.
    check_gene_link_to_mrna(\%omniscient);

    #check level1 has subfeature else we remove it
  	remove_orphan_l1(\%omniscient);

    # To keep track of How many Warnings we got....
    foreach my $thematic (keys %WARNS){
  		my $nbW = $WARNS{$thematic};
  		if($nbW > $nbWarnLimit){
  			print "Actually we had $nbW warning message: $thematic\n";
  		}	
  	}

  	print "Parsing done in ", time() - $start_run," seconds\n\n\n";
  	
    #return
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

# The method read a gff3 feature (one line), Check for correctly formated gff3 feature (ID and Parent attribute) 
sub manage_one_feature{
	
	my ($feature, $omniscient, $mRNAGeneLink, $duplicate, $miscCount)=@_;

		my $seq_id = $feature->seq_id;					#col1
		my $source_tag = lc($feature->source_tag);		#col2
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
	    if ($primary_tag eq "gene" or $primary_tag eq "sts" ) {
	    	#get ID
	    	if($feature->has_tag('ID')){
	    		$id = lc($feature->_tag_value('ID')) ;
	    	}
	    	else{warn "gff3 reader error level2: No ID attribute found @ for the feature: ".$feature->gff_string()."\n";
	    		$id = _manage_ID($omniscient, $duplicate, $miscCount, $feature, 'level1', $primary_tag);
			}
    	
    		#Save feature, but first check for duplicated line/feature and duplicated ID.
    		if(! exists_keys($omniscient,('level1',$primary_tag,$id))){
	        	$omniscient->{"level1"}{$primary_tag}{$id}=$feature;
	        }
	        else{ # ID is duplicated
	        	if(_feature_is_duplicated( $omniscient->{'level1'}{$primary_tag}{$id}, $feature, $miscCount, $primary_tag) ) {	
	        	# line/feature is duplicated
		        	push (@{$duplicate->{"level1"}{$primary_tag}{$id}}, $feature);
				}
	        }
	    }

      	###################################################
      	########## Manage feature WITHout CHILD  ##########		== LEVEL3 ==
      	###################################################

      	elsif ( ($primary_tag eq "cds") or ($primary_tag eq "exon") or ($primary_tag eq "stop_codon") or ($primary_tag eq "start_codon") or ($primary_tag eq "three_prime_utr") or ($primary_tag eq "five_prime_utr") or ($primary_tag eq "utr"),
      	 or ($primary_tag eq "selenocysteine") or ($primary_tag eq "non_canonical_three_prime_splice_site") or ($primary_tag eq "non_canonical_five_prime_splice_site") or ($primary_tag eq "stop_codon_read_through"),
      	 or ($primary_tag eq "sig_peptide") or ($primary_tag eq "tss") or ($primary_tag eq "tts") or ($primary_tag eq "intron") ){

      		# manage ID
      		if(! $feature->has_tag('ID')){
      			warn "WARNING gff3 reader level3: No ID attribute found @ for the feature: ".$feature->gff_string()."\n";
				$id = _manage_ID($omniscient, $duplicate, $miscCount, $feature, 'level3', $primary_tag);	
			}

			# manage Parent
			my @parentList;
      		if($feature->has_tag('Parent')){
      			@parentList = $feature->get_tag_values('Parent');
			}
			else{ # In that case we create a uniq parent ID to create a preper omniscient structure. But the feature itself stay intact without parentID.
				warn "WARNING gff3 reader level3: No Parent attribute found @ for the feature: ".$feature->gff_string()."\n";
				$miscCount->{'noParent'}{$primary_tag}++;
				$parent = $primary_tag."-".$miscCount->{'noParent'}{$primary_tag};
				create_or_replace_tag($feature,'Parent',$parent); #add parent attribute
				push(@parentList, $parent);
			}
			
			#Save feature and check duplicates	(treat also cases where there is multiple parent. => In that case we expand to create a uniq feature for each)
			my $cptParent=0; # to check if it is a multiple parent case
      		my $allParent = scalar @parentList;
      		foreach my $parent (@parentList){ # first feature level3 with this primary_tag linked to the level2 feature
				$cptParent++;

				#Level3 key doesn't exist
				if(! exists_keys($omniscient,('level3',$primary_tag,lc($parent)))){
					
					# It is a multiple parent case
					if($allParent > 1){
						# Not the first parent, we have to clone the feature !!
						if($cptParent > 1){ 
							my $feature_clone=clone($feature);
							create_or_replace_tag($feature_clone,'Parent',$parent); #modify Parent To keep only one

							#As we cloned the feature, we haev to take care of ID to have uniq one (case of exon,etc...)
							if( ($primary_tag ne "cds") and (index($primary_tag, 'utr') == -1) ){ #avoid case where we can have the same ID (multi-feature)
								my $substr = $feature->_tag_value('Parent');
								my $clone_id=$feature->_tag_value('ID');
								$clone_id=~ s/$substr//;
								if($clone_id eq $feature->_tag_value('ID')){#Substring didn't work
									$clone_id=$feature->_tag_value('ID')."$cptParent";
								}else{$clone_id="$parent$clone_id";}
								create_or_replace_tag($feature_clone,'ID',$clone_id); #modify Parent To keep only one
							}

							push (@{$omniscient->{"level3"}{$primary_tag}{lc($parent)}}, $feature_clone);
						}
						# It is the first parent. Do not clone the feature
						else{
							create_or_replace_tag($feature,'Parent',$parent); #modify Parent To keep only one
							push (@{$omniscient->{"level3"}{$primary_tag}{lc($parent)}}, $feature);
						}
					}
					else{ #the simpliest case. One parent only
						push (@{$omniscient->{"level3"}{$primary_tag}{lc($parent)}}, $feature);
					}
				}
				#Level3 key exists
				else{  # If not the first feature level3 with this primary_tag linked to the level2 feature
					# check among list of feature level3 already exits with an identical ID.
					my $is_dupli=undef;
					foreach my $feat ( @{$omniscient->{"level3"}{$primary_tag}{lc($parent)} }){
						# case where same feature is spread on different location (i.e utr, cds). In that case, following the gff3 gene ontologie specification, the ID can be share by all "peace of feature" building the feature entity.
						if(($primary_tag eq "cds") or ($primary_tag =~ /utr/)){
							if(lc($feat->seq_id().$feat->start().$feat->end().$feat->_tag_value('ID')) eq  lc($feature->seq_id().$feature->start().$feature->end().$feature->_tag_value('ID'))){
								my $id = $feature->seq_id().$feature->start().$feature->end().$feature->_tag_value('ID');
								push ( @{$duplicate->{"level3"}{$primary_tag}{lc($id )}}, $feature );
								$is_dupli=1;
								last;
							}
						}# case where a feature are an uniq element (exon, etc.)
						elsif(lc($feat->_tag_value('ID')) eq lc($feature->_tag_value('ID'))){
							my $id = $feature->_tag_value('ID');
							push ( @{$duplicate->{"level3"}{$primary_tag}{lc($id)}}, $feature );
							$is_dupli=1;
							last;
						}
					}
					#It is not a duplicated feature => save it in omniscient
					if(! $is_dupli){
						# It is a multiple parent case
						if($allParent > 1){
							# Not the first parent, we have to clone the feature !!
							if($cptParent > 1){ 
								my $feature_clone=clone($feature);
								create_or_replace_tag($feature_clone,'Parent',$parent); #modify Parent To keep only one
								
								#Take care of ID to have uniq one (case of exon,etc...)
								if( ($primary_tag ne "cds") and (index($primary_tag, 'utr') == -1) ){
									my $substr = $feature->_tag_value('Parent');
									my $clone_id=$feature->_tag_value('ID');
									$clone_id=~ s/$substr//;
									if($clone_id eq $feature->_tag_value('ID')){#Substring didn't work
										$clone_id=$feature->_tag_value('ID')."$cptParent";
									}else{$clone_id="$parent$clone_id";}
									create_or_replace_tag($feature_clone,'ID',$clone_id); #modify Parent To keep only one
								}

								push (@{$omniscient->{"level3"}{$primary_tag}{lc($parent)}}, $feature_clone);
							}
							# It is the first parent. Do not clone the feature
							else{
								create_or_replace_tag($feature,'Parent',$parent); #modify Parent To keep only one
								push (@{$omniscient->{"level3"}{$primary_tag}{lc($parent)}}, $feature);
							}
						}
						else{ #the simpliest case. One parent only
							push (@{$omniscient->{"level3"}{$primary_tag}{lc($parent)}}, $feature);
						}
					}
				}
      		}
      	}

      	##############################################
      	########## Manage feature the rest  ##########		== LEVEL2 ==
      	##############################################
      	elsif( ($primary_tag eq "mrna") or ($primary_tag eq "ncrna") or ($primary_tag eq "mirna") or ($primary_tag eq "lcrna") or ($primary_tag eq "rrna") or ($primary_tag eq "srp_rna"),
      		or ($primary_tag eq "snrna") or ($primary_tag eq "lincrna") or ($primary_tag eq "trna") or ($primary_tag eq "trna_pseudogene") or ($primary_tag eq "snorna") or ($primary_tag eq "misc_rna"),
      		or ($primary_tag eq "rnase_p_rna") or ($primary_tag eq "tmrna")){
    		
    		#get ID
	    	if($feature->has_tag('ID')){
		    	$id = lc($feature->_tag_value('ID'));
			}

			else{warn "gff3 reader error level2: No ID attribute found @ for the feature: ".$feature->gff_string()."\n";
				$id = _manage_ID($omniscient, $duplicate, $miscCount, $feature, 'level2', $primary_tag);
			}

			#get Parent
			if($feature->has_tag('Parent')){
				$parent = lc($feature->_tag_value('Parent'));
			}
			else{warn "WARNING gff3 reader level2 : No Parent attribute found for @ the feature: ".$feature->gff_string()."\n";
				$miscCount->{'noParent'}{$primary_tag}++;
				$parent = $primary_tag."-".$miscCount->{'noParent'}{$primary_tag};	
				$feature->add_tag_value('Parent', $parent);
			}

			# keep track of link between level2->leve1
  			if (! exists ($mRNAGeneLink->{lc($id)})){ 
				$mRNAGeneLink->{lc($id)}=lc($parent);
	 		}

	 		#Save feature and check duplicates
	 		if(! exists_keys($omniscient,('level2',$primary_tag,lc($parent)))){# case of first feature l2 linked to the level1 feature 
      			push (@{$omniscient->{"level2"}{$primary_tag}{lc($parent)}}, $feature);
      		}
      		else{# case where isoforms exist

      			# check among list of feature if one with a similar ID exists.
      			if(_feature_is_duplicated(\@{$omniscient->{"level2"}{$primary_tag}{lc($parent)}}, $feature, $miscCount, $primary_tag)){
					push (@{$duplicate->{"level2"}{$primary_tag}{lc($parent)}}, $feature);
     			}
     			else{ # No duplicate found, we keep it as isoform
	     			push (@{$omniscient->{"level2"}{$primary_tag}{lc($parent)}}, $feature);
	     		}		
      		}
      	}

      	#####################################
	    ########## Manage REPEATS // Protein2genome // stuff in match/matchpart ###########	/!\ NO level3 features /!\ Compare to gene stuff we replace primary_tag by source_tag when filling the omniscient hash
	    #####################################
      	elsif($primary_tag =~ /match/) { #If primary tag contain match we are in case where we should take in consideration $source_tag instead $primary_tag

      		
      		#									== LEVEL1 == 
			if ($primary_tag ne "match_part" and $primary_tag ne "similarity") { #means we are in level1 case (From maker it could be match, protein_match)
	      		#get ID
		    	if($feature->has_tag('ID')){
		    		$id = lc($feature->_tag_value('ID')) ;
		    	}
		    	else{warn "gff3 reader error level2: No ID attribute found @ for the feature: ".$feature->gff_string()."\n";
					$id = _manage_ID($omniscient, $duplicate, $miscCount, $feature, 'level1', $source_tag);
				}
	    	
	    		#Save feature
	    		if(! exists_keys($omniscient,('level1',$source_tag,$id))){
		        	$omniscient->{"level1"}{$source_tag}{$id}=$feature;
		        }
		        else{push (@{$duplicate->{"level1"}{$source_tag}{$id}}, $feature);}
	    	}
	    	
	    	      #								== LEVEL2 == 
			else { # means we are in a case where ($primary_tag eq "match_part" or $primary_tag eq "similarity") => and they are level2

				#get ID
		    	if($feature->has_tag('ID')){
			    	$id = lc( $feature->get_tag_values('ID') );
				}
				else{warn "gff3 reader error level2: No ID attribute found @ for the feature: ".$feature->gff_string()."\n";
					$id = _manage_ID($omniscient, $duplicate, $miscCount, $feature, 'level2', $source_tag);
				}

				#get Parent
				if($feature->has_tag('Parent')){
					$parent = lc($feature->_tag_value('Parent'));
				}
				else{warn "WARNING gff3 reader level2 : No Parent attribute found for @ the feature: ".$feature->gff_string()."\n";
					$miscCount->{'noParent'}{$source_tag}++;
					$parent = $source_tag."-".$miscCount->{'noParent'}{$source_tag};
					$feature->add_tag_value('Parent', $parent);	
				}


	  			if (! exists ($mRNAGeneLink->{$id})){ # keep track of link between level2->leve1
					$mRNAGeneLink->{$id}=$parent;
		 		}

		 		#Save feature and check duplicates
		 		if(! exists_keys($omniscient,('level2',$source_tag,lc($parent)))){ # case of first feature l2 linked to the level1 feature 
	      			push (@{$omniscient->{"level2"}{$source_tag}{lc($parent)}}, $feature);
	      		}
	      		else{ # case where isoforms exist
	      			# check among list of feature if one with a similar ID exists.
	      			if(_feature_is_duplicated(\@{$omniscient->{"level2"}{$source_tag}{lc($parent)}}, $feature, $miscCount, $source_tag)){
	      				push (@{$duplicate->{"level2"}{$source_tag}{lc($parent)}}, $feature);
	      			}
	      			else{ # No similar ID found, we keep it as isoform
	     				push (@{$omniscient->{"level2"}{$source_tag}{lc($parent)}}, $feature);
	     			}		
	      		}
			}
      	}
      	else{
      		print "gff3 reader warning: $primary_tag still not taken in account ! Please modify the code to define on of the three level of this feature.\n";
      	}	

}

##==============================================================

# This method is called when no ID exist for the feature. It creates one with a prefix and number and add it to the attribute ID.
# This check if the feature already exists without ID but with same parent (when parent exists)
sub _manage_ID{
	my	($omniscient, $duplicate, $miscCount, $feature, $level, $tag)=@_;

	#Before to create the ID we check if line is a duplicate Else we just remove it.
	#####
	#get Parent
	my $parent;
	if($feature->has_tag('Parent')){
		$parent = lc($feature->_tag_value('Parent'));
	}
	else{
		if($level eq 'level1'){
			$parent = "undef"; 
		}
		else{
			warn "WARNING _manage_ID : No Parent attribute found for @ the feature: ".$feature->gff_string()."\n";
		}
	}
	if($parent){
		#Check if feature duplicated
		if(exists_keys($omniscient,($level,$tag,lc($parent)))){
			if(_feature_is_duplicated(\@{$omniscient->{$level}{$tag}{lc($parent)}}, $feature, $miscCount, $tag)){
				#feature identic already exists
			     push (@{$duplicate->{$level}{$tag}{lc($parent)}}, $feature);
			}
			else{ # Similar ID found,  No similar feature found
			    my $feature_clone=clone($feature);
			    _create_ID($miscCount, $feature, $tag);
			    warn "WARNING _manage_ID : Duplicate ID found @ for the feature: ".$feature_clone->gff_string().". We change it !\n";
			}
		}
		# No similar ID found
		else{ #create the ID
			_create_ID($miscCount, $feature, $tag);
		}
	}
	else{ #create the ID - Case where no parent
		_create_ID($miscCount, $feature, $tag);
	}
}

# Check if a feature is a duplicate
# To do that check we compare position (start, end and chr/scaffold/contig) and the ID.
# If ID is absent we compare the Parent attribute instead, or none of them if both are absent.
# $list is a feature of a reference of list of feature.
sub _feature_is_duplicated{
	my ($list, $feature, $miscCount, $tag)=@_;

	#check if list is a feature or a reference of a list of feature
	my @list_feature; # will be a list of feature.
	if(ref($list) eq 'ARRAY'){
		@list_feature=@$list; #it's a reference of a list of feature
	}
	else{push (@list_feature, $list);} # it's a feature

	my $is_dupli=undef;

	#current data for testing duplicate
	my $current_string_line_dupli;
	if(! $feature->has_tag('ID')){
		if($feature->has_tag('Parent')){
			$current_string_line_dupli = $feature->seq_id().$feature->start().$feature->end().$feature->_tag_value('Parent');
		}
		else{
			$current_string_line_dupli = $feature->seq_id().$feature->start().$feature->end(); #case level1
		}
	}
	else{
		$current_string_line_dupli = $feature->seq_id().$feature->start().$feature->end().$feature->_tag_value('ID'); 
	}

	#Check all the level2 list element
	foreach my $feature_in_omniscient ( @list_feature){

		#get data in omniscient
		my $string_line_dupli_in_omniscient;
		if(! $feature->has_tag('ID')){
			if($feature_in_omniscient->has_tag('Parent')){
				$string_line_dupli_in_omniscient = $feature_in_omniscient->seq_id().$feature_in_omniscient->start().$feature_in_omniscient->end().$feature_in_omniscient->_tag_value('Parent');
			}
			else{
				$string_line_dupli_in_omniscient = $feature_in_omniscient->seq_id().$feature_in_omniscient->start().$feature_in_omniscient->end();
			}
		}
		else{
		 	$string_line_dupli_in_omniscient = $feature_in_omniscient->seq_id().$feature_in_omniscient->start().$feature_in_omniscient->end().$feature_in_omniscient->_tag_value('ID');
		}
		
		#compare data if identical feature

		if($current_string_line_dupli eq $string_line_dupli_in_omniscient){
			$is_dupli=1;
			last;
		}
		#compare data if identical ID. Same ID but not same feature, so we just change the ID
		if($feature_in_omniscient->has_tag('ID') and $feature->has_tag('ID')){
			if($feature_in_omniscient->_tag_value('ID') eq $feature->_tag_value('ID')){
				warn "WARNING_T _manage_ID : Duplicate ID found @ for the feature ".$feature->gff_string().". We change it !\n";
				_create_ID($miscCount, $feature, $tag)
			}
		}
	}
	return $is_dupli;
}

# create the ID and add it to the feature.
sub _create_ID{
	my	($miscCount, $feature, $tag)=@_;
	
	$miscCount->{'noID'}{$tag}++;
	my $id = $tag."-".$miscCount->{'noID'}{$tag};
	create_or_replace_tag($feature,'ID', $id);
	
	return $id;
}

# check if mrNA have is PArenttal gene existing. If not we create it.
sub check_gene_link_to_mrna{
	my ($hash_omniscient)=@_;

	foreach my $primary_tag_l2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
		foreach my $id_l1 (keys %{$hash_omniscient->{'level2'}{$primary_tag_l2}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
			my $l1_exist=undef;
			foreach my $primary_tag_l1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
				if(exists ($hash_omniscient->{'level1'}{$primary_tag_l1}{$id_l1})){
					$l1_exist=1;
					last;
				}
			}
			if(! $l1_exist){
				warn "WARNING gff3 reader level2 : No Parent feature found with the ID @ ".$id_l1.". We will create one.\n";
				my $gene_feature=clone($hash_omniscient->{'level2'}{$primary_tag_l2}{$id_l1}[0]);#create a copy of the first mRNA feature;
				my $new_ID = $gene_feature->_tag_value('Parent');
				create_or_replace_tag($gene_feature,'ID', $new_ID); #modify ID to replace by parent value
				$gene_feature->remove_tag('Parent'); # remove parent ID because, none.
				check_level1_positions($hash_omniscient, $gene_feature);	# check start stop if isoforms exists
				
				#Deal case where we reconstruct other thing than a gene
				my $primary_tag_l1=undef;
				if(lc($gene_feature->primary_tag) =~ /match/){ $primary_tag_l1="match"; }
				else{ $primary_tag_l1="gene"; }

				$gene_feature->primary_tag($primary_tag_l1); # change primary tag
				$hash_omniscient->{"level1"}{$primary_tag_l1}{lc($new_ID)}=$gene_feature; # now save it in omniscient
				#print "feature level1 created: ".$gene_feature->gff_string."\n";
			}
		}
	}
}


#
sub modelate_utr_and_cds_features_from_exon_features_and_cds_start_stop{

	my ($exon_features, $ORFstart, $ORFend)=@_;

	my @cds_features;
	my @utr3_features;
	my @utr5_features;
	my $strand = $exon_features->[0]->strand;
  	my @exon_features_sorted = sort {$a->start <=> $b->start} @{$exon_features}; # be sure that exon list is sorted   

  	my $cds_counter=1;
  	my $utr3_counter=1;
  	my $utr5_counter=1;
 	foreach my $exon_feature (@exon_features_sorted){

	    # exon overlap fully a CDS
	    if( ($exon_feature->end >= $ORFend) and ($exon_feature->start <= $ORFstart) ){

 			my $cds_feature=clone($exon_feature);#create a copy of the feature 					exon    ====================================
 			$cds_feature->start($ORFstart); #modify start 											 cds     ============================ 
 			$cds_feature->end($ORFend); #modify end
 			$cds_feature->primary_tag('CDS');
 			#get old name
 			my $ID = $cds_feature->_tag_value('ID');
 			create_or_replace_tag($cds_feature,'ID',$ID.'-cds-'.$cds_counter); #modify name
 			push(@cds_features, $cds_feature);#save that cds
 			$cds_counter++;
 			if($exon_feature->start < $ORFstart){
 				my $utr_feature=clone($exon_feature);#create a copy of the feature 
 				$utr_feature->end($ORFstart-1); #modify start 
 				if ( ($strand == -1) or ($strand eq "-") ) {
 					$utr_feature->primary_tag('three_prime_UTR');
 					create_or_replace_tag($utr_feature,'ID',$ID.'-utr3-'.$utr3_counter); #modify name
	 				push(@utr3_features, $utr_feature);#save that cds
	 				$utr3_counter++;
	 			}else{
	 				$utr_feature->primary_tag('five_prime_UTR');
	 				create_or_replace_tag($utr_feature,'ID',$ID.'-utr5-'.$utr5_counter); #modify name
	 				push(@utr5_features, $utr_feature);#save that cds
	 				$utr5_counter++;
	 			}
 			}
 			if($exon_feature->end > $ORFend){
 				my $utr_feature=clone($exon_feature);#create a copy of the feature 
 				$utr_feature->start($ORFend+1); #modify start 
 				if ( ($strand == -1) or ($strand eq "-") ) {
 					$utr_feature->primary_tag('five_prime_UTR');
	 				create_or_replace_tag($utr_feature,'ID',$ID.'-utr5-'.$utr5_counter); #modify name
	 				push(@utr5_features, $utr_feature);#save that cds
	 				$utr5_counter++;
	 			}else{
	 				$utr_feature->primary_tag('three_prime_UTR');
 					create_or_replace_tag($utr_feature,'ID',$ID.'-utr3-'.$utr3_counter); #modify name
	 				push(@utr3_features, $utr_feature);#save that cds
	 				$utr3_counter++;
	 			}
 			}
		}
		# cds overlap fully an exon
		elsif( ($exon_feature->end <= $ORFend) and ($exon_feature->start >= $ORFstart) ){
 			my $cds_feature=clone($exon_feature);#create a copy of the feature 						exon    ========================
 			$cds_feature->primary_tag('CDS');
 			#get old name 																			cds  ===============================
 			my $ID = $cds_feature->_tag_value('ID');
 			create_or_replace_tag($cds_feature,'ID',$ID.'-cds-'.$cds_counter); #modify name
 			push(@cds_features, $cds_feature);#save that cds
 			$cds_counter++;
		}
		# cds overp partially an exon
	    elsif( ($exon_feature->end >= $ORFstart) and ($exon_feature->start <= $ORFend) ){ #they overlap
	      
	      if($exon_feature->start >= $ORFstart){ # cds overlap start of exon                                    exon ===============================
	      	#Manage CDS
	      	my $cds_feature=clone($exon_feature);#create a copy of the feature 						cds ===============================
 			$cds_feature->end($ORFend); #modify end
 			$cds_feature->primary_tag('CDS');
 				#get old name
 			my $ID = $cds_feature->_tag_value('ID');
 			create_or_replace_tag($cds_feature,'ID',$ID.'-cds-'.$cds_counter); #modify name
 			push(@cds_features, $cds_feature);#save that cds
 			$cds_counter++;
 			#manage UTR
 			my $utr_feature=clone($exon_feature);#create a copy of the feature 
 			$utr_feature->start($ORFend+1); #modify end 
 			$ID = $utr_feature->_tag_value('ID');
	 		if ( ($strand == -1) or ($strand eq "-") ) {
	 			$utr_feature->primary_tag('five_prime_UTR');
	 			create_or_replace_tag($utr_feature,'ID',$ID.'-utr5-'.$utr5_counter); #modify name
	 			push(@utr5_features, $utr_feature);#save that cds
	 			$utr5_counter++;

	 		}else{
				$utr_feature->primary_tag('three_prime_UTR');
	 			create_or_replace_tag($utr_feature,'ID',$ID.'-utr3-'.$utr3_counter); #modify name
	 			push(@utr3_features, $utr_feature);#save that cds
	 			$utr3_counter++;
	 		}

	      }
	      else{ #cds overlap start end exon
	       	#Manage CDS
	       	my $cds_feature=clone($exon_feature);#create a copy of the feature
 			$cds_feature->start($ORFstart); #modify start 										exon ===============================
 			$cds_feature->primary_tag('CDS');	
 				#get old name 																					cds =====================================
 			my $ID = $cds_feature->_tag_value('ID');
 			create_or_replace_tag($cds_feature,'ID',$ID.'-cds-'.$cds_counter); #modify name
 			push(@cds_features, $cds_feature);#save that cds
 			$cds_counter++;
	 		 #Manage UTR
 			my $utr_feature=clone($exon_feature);#create a copy of the feature 
 			$utr_feature->end($ORFstart-1); #modify start 
 			$ID = $utr_feature->_tag_value('ID');
	 		if ( ($strand == -1) or ($strand eq "-") ) {
	 			$utr_feature->primary_tag('three_prime_UTR');
	 			create_or_replace_tag($utr_feature,'ID',$ID.'-utr3-'.$utr3_counter); #modify name
	 			push(@utr3_features, $utr_feature);#save that cds
	 			$utr3_counter++;
	 		}else{
	 			$utr_feature->primary_tag('five_prime_UTR');
	 			create_or_replace_tag($utr_feature,'ID',$ID.'-utr5-'.$utr5_counter); #modify name
	 			push(@utr5_features, $utr_feature);#save that cds
	 			$utr5_counter++;
	 		}
	      }
	    }###### Only UTR part
	    else{ #Does not overlap
	    	if($exon_feature->end < $ORFstart){ #UTR5 in + strand
		    	my $utr_feature=clone($exon_feature);#create a copy of the feature 			exon ===============================
	 			#get old name 																											 cds ===============================
	 			my $ID = $utr_feature->_tag_value('ID');
	 			if ( ($strand == -1) or ($strand eq "-") ) {
	 				$utr_feature->primary_tag('three_prime_UTR');
	 				create_or_replace_tag($utr_feature,'ID',$ID.'-utr3-'.$utr3_counter); #modify name
	 				push(@utr3_features, $utr_feature);#save that cds
	 				$utr3_counter++;
	 			}else{
	 				$utr_feature->primary_tag('five_prime_UTR');
	 				create_or_replace_tag($utr_feature,'ID',$ID.'-utr5-'.$utr5_counter); #modify name
	 				push(@utr5_features, $utr_feature);#save that cds
	 				$utr5_counter++;
	 			}
	 			
	 		}
	 		else{	#UTR3 in + strand
		    	my $utr_feature=clone($exon_feature);#create a copy of the feature 													exon ===============================
	 			#get old name
	 			my $ID = $utr_feature->_tag_value('ID'); 									#cds ===============================
	 			if ( ($strand == -1) or ($strand eq "-") ) {
	 				$utr_feature->primary_tag('five_prime_UTR');
	 				create_or_replace_tag($utr_feature,'ID',$ID.'-utr5-'.$utr5_counter); #modify name
	 				push(@utr5_features, $utr_feature);#save that cds
	 				$utr5_counter++;
	 			}else{
	 				$utr_feature->primary_tag('three_prime_UTR');
	 				create_or_replace_tag($utr_feature,'ID',$ID.'-utr3-'.$utr3_counter); #modify name
	 				push(@utr3_features, $utr_feature);#save that cds
	 				$utr3_counter++;
	 			}	
	 		}
	    }
 	}
my @utr5_features_sorted=sort {$a->start <=> $b->start} @utr5_features;
my @cds_features_sorted=sort {$a->start <=> $b->start} @cds_features;
my @utr3_features_sorted=sort {$a->start <=> $b->start} @utr3_features;
return \@utr5_features_sorted, \@cds_features_sorted, \@utr3_features_sorted; #really utr5 and utr3 that are return
}


# Check the start and end of mRNA based a list of feature like list of exon;
sub check_mrna_positions{

  my ($mRNA_feature, $exon_list)=@_;

  my @exon_list_sorted = sort {$a->start <=> $b->start} @{$exon_list};
  ######
  #Modify mRNA start-end based on exon features
  my $exonStart=$exon_list_sorted[0]->start;
  my $exonEnd=$exon_list_sorted[$#exon_list_sorted]->end;
  #check start
  if ($mRNA_feature->start != $exonStart){
    $mRNA_feature->start($exonStart);
  }
  #check stop
  if($mRNA_feature->end != $exonEnd){
    $mRNA_feature->end($exonEnd);
  }
}

# @Purpose: Remove the level1 feature that havn't subfeature linked to it
# @input: 1 => hash(omniscient hash)
# @output: none
sub remove_orphan_l1{
	my ($hash_omniscient)=@_;

	foreach my $tag_l1 (keys %{$hash_omniscient->{'level1'}}){
	  	foreach my $id_l1 (keys %{$hash_omniscient->{'level1'}{$tag_l1}}){
	     
		    my $neverfound="yes";
		    foreach my $tag_l2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
		        if ( exists ($hash_omniscient->{'level2'}{$tag_l2}{$id_l1} ) ){
		          $neverfound=undef;last
		        }   
		    }
		    if($neverfound){
		    	warn "WARNING gff3 reader level1 : No child feature found for the $tag_l1 with ID @ ".$id_l1.". We remove it.\n";
		    	delete $hash_omniscient->{'level1'}{$tag_l1}{$id_l1}; # delete level1
		    }
	 	}
	}
}

1;
