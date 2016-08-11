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

use constant LEVEL1 => { "gene" => 1, "sts" => 2, "match" => 3 };
use constant LEVEL2 => { "mrna" => 1, "ncrna" => 2, "mirna" => 3, "lcrna" => 4, "rrna" => 5, "srp_rna" => 6, "snrna" => 7, "lincrna" => 8, "trna" => 9, "trna_pseudogene" => 10, "snorna" => 11, "misc_rna" => 12, "rnase_p_rna" => 13, "tmrna" => 14, "match_part" => 15, "similarity" => 16, "rna" => 17};
use constant LEVEL3 => { "cds" => 1, "exon" => 2, "stop_codon" => 3, "start_codon" => 4, "three_prime_utr" => 5, "five_prime_utr" => 6, "utr" => 7, "selenocysteine" => 8, "non_canonical_three_prime_splice_site" => 8, "non_canonical_five_prime_splice_site" => 10,
						"stop_codon_read_through" => 11, "sig_peptide" => 12, "tss" => 13, "tts" => 14, "intron" => 15 };
use constant SPREADFEATURE => {"cds" => 1, "three_prime_utr" => 2, "five_prime_utr" => 3, "utr" => 4};

# Save in omniscient hash (sorted in a specific way (3 levels)) a whole gff3 file
sub slurp_gff3_file_JD {
	
	my ($self, $file) = @_  ;
	
	my $gffio = Bio::Tools::GFF->new(-file => $file, -gff_version => 3);	


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
	my %uniqID;# Hash to follow up with an uniq identifier every feature
	my %infoSequential;# Hash to store any counter. Will be use to create a new uniq ID
	my $comonTag=undef;
	my $last_l1_f=undef;
	my $last_l2_f=undef;

	#read every lines
	while( my $feature = $gffio->next_feature()) {
		($comonTag, $last_l1_f, $last_l2_f) = manage_one_feature($feature, \%omniscient, \%mRNAGeneLink, \%duplicate, \%miscCount, \%uniqID, \%infoSequential, $comonTag, $last_l1_f, $last_l2_f);
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

    #check level1 has subfeature else we remove it
  	remove_orphan_l1(\%omniscient, \%miscCount, \%uniqID); #or fix if level2 is missing (refseq case)

#   #Check sequential if we can fix cases
#   check_sequential(\%infoSequential, \%omniscient, \%miscCount);

#    #Check relationship between mRNA and gene.
#    check_gene_link_to_mrna(\%omniscient);



    # To keep track of How many Warnings we got....
    foreach my $thematic (keys %WARNS){
  		my $nbW = $WARNS{$thematic};
  		if($nbW > $nbWarnLimit){
  			print "Actually we had $nbW warning message: $thematic\n";
  		}	
  	}

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
	
	my ($feature, $omniscient, $mRNAGeneLink, $duplicate, $miscCount, $uniqID, $infoSequential, $comonTag, $last_l1_f, $last_l2_f)=@_;

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
		my @comonTagList=('locus_tag'); #Tag used in old gff format allowing to group features together. Priority to comonTag compare to sequential read 
		my $currentBucket=undef;

		####################################################
	    ########## Manage feature WITHOUT parent ###########	== LEVEL1 ==
	    ####################################################
	    if (exists(LEVEL1->{$primary_tag}) ) {
	    	
	    	################
	    	# REPEATS case #
	    	if($primary_tag eq "match"){ #Manage REPEATS // Protein2genome // stuff in match/matchpart ###########	/!\ NO level3 features /!\ Compare to gene stuff we replace primary_tag by source_tag when filling the omniscient hash
	    		$primary_tag=$source_tag;
	    	}

	    	##########
			# get ID #
	    	$id = lc(_check_uniq_id($miscCount, $uniqID, $feature));

	    	#####################
	    	# Ckeck duplication #
    		if(! _it_is_duplication($duplicate, $omniscient, $uniqID, $feature)){
    			
	    		################
	    		# Save feature #
	    		$last_l1_f = $feature;
    			$omniscient->{"level1"}{$primary_tag}{$id}=$feature;
 	       		
 	       		#############
 	       		# COMON TAG #
		        my $comonTag =_comon_tag(\@comonTagList, $feature);

				if($comonTag){
		    		$infoSequential->{$comonTag}{'level1'}=$id;
		    		$last_l1_f=undef; #opposite to I have a comon tag
		    		$last_l2_f=undef; #opposite to I have a comon tag
		    	}
		    	else{
		    		$comonTag=undef; #reinitialize to empty.
		    	}
		    }
	    }

      	###################################################
      	########## Manage feature WITHout CHILD  ##########		== LEVEL3 ==
      	###################################################
      	elsif ( exists(LEVEL3->{$primary_tag}) ){

      		##########
			# get ID #
	    	$id = lc(_check_uniq_id($miscCount, $uniqID, $feature));

	    	##############
			# get Parent #
			my @parentList;
      		if($feature->has_tag('Parent')){
      			@parentList = $feature->get_tag_values('Parent');
			}
			else{ # In that case we create a uniq parentID to create a proper omniscient structure. But the feature itself stay intact without parentID.
				warn "WARNING gff3 reader level3: No Parent attribute found @ for the feature: ".$feature->gff_string()."\n";
				
				######################
				# NEED THE LEVEL2 ID #
				my $l2_id="";
				# case where No level2 feature defined yet - I will need a bucketL2
				if(! $last_l2_f ){ 
					$l2_id = _create_ID($miscCount, $uniqID, $feature, $id, "nbis_noL2id");					
					$last_l2_f = clone($feature);
					create_or_replace_tag($last_l2_f,'ID',$l2_id); #modify Parent To keep only one
					$last_l1_f->primary_tag('rna');
				}
				else{ # case where previous level2 exists
					$l2_id=$last_l2_f->_tag_value('ID');
				}

				#############
 	       		# COMON TAG #
				my $comonTag =_comon_tag(\@comonTagList, $feature);

				if($comonTag){ #Previous Level up feature had a comon tag
			
		    		push( @{$infoSequential->{$comonTag}{$l2_id}{'level3'}}, $feature );			
			    	return $comonTag, $last_l1_f, $last_l2_f;								#### STOP HERE AND RETURN
				}
				else{# No comon tag found 
					######################
					# NEED THE LEVEL1 ID #
					my $l1_id="";			
					if($last_l1_f){ # case where previous level1 exists
						$l1_id=$last_l1_f->_tag_value('ID');
					}
					else{ # case where No level2 feature defined yet - I will need a bucketL2
						$l1_id = _create_ID($miscCount, $uniqID, $feature, $id, "nbis_noL1id");					
						$last_l1_f = clone($feature);
						create_or_replace_tag($last_l1_f,'ID',$l1_id); #modify Parent To keep only one
						$last_l1_f->primary_tag('gene');
					}

					push( @{$infoSequential->{$$l1_id}{$l2_id}{'level3'}}, $feature );
					return $l2_id, $last_l1_f, $last_l2_f;								#### STOP HERE AND RETURN
				}
			}
			
			####################
			# HANDLE PARENT(S) # #Save feature and check duplicates	(treat also cases where there is multiple parent. => In that case we expand to create a uniq feature for each)
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

							#As we cloned the feature, we have to take care of ID to have uniq one (case of exon,etc...)
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
					foreach my $feat ( @{ $omniscient->{"level3"}{$primary_tag}{lc($parent)} } ) {
						# case where same feature is spread on different location (i.e utr, cds). In that case, following the gff3 gene ontologie specification, the ID can be share by all "peace of feature" building the feature entity.
						if( _it_is_duplication($duplicate, $omniscient, $uniqID, $feature) ){
								$is_dupli=1;
								last;
						}# case where a feature are an uniq element (exon, etc.)
						elsif( lc($feat->_tag_value('ID')) eq lc($feature->_tag_value('ID')) ){
							my $id = $feature->_tag_value('ID');

							if( $feat->start() == $feature->start() and $feat->end() == $feature->end() ){
								push ( @{$duplicate->{"level3"}{$primary_tag}{lc($id)}}, $feature );
								$is_dupli=1;
								last;
							}
							else{ #case we have to change name ! GFF3 NON CONFORME
								$miscCount->{$id}++;
								$id=$id."-".$miscCount->{$id};
								create_or_replace_tag($feature,'ID',$id); #modify Parent To keep only one
							}
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
      	elsif ( exists(LEVEL2->{$primary_tag}) ) {
    		if ($primary_tag eq "match_part" or $primary_tag eq "similarity") { ########## Manage REPEATS // Protein2genome // stuff in match/matchpart ###########	/!\ NO level3 features /!\ Compare to gene stuff we replace primary_tag by source_tag when filling the omniscient hash       	
				$primary_tag=$source_tag;  
    		}

    		##########
			# get ID #
	    	$id = lc(_check_uniq_id($miscCount, $uniqID, $feature));

			##############
			# get Parent #
			if($feature->has_tag('Parent')){
				$parent = lc($feature->_tag_value('Parent'));
			}
			else{warn "WARNING gff3 reader level2 : No Parent attribute found for @ the feature: ".$feature->gff_string()."\n";

				#############
 	       		# COMON TAG #
				my $comonTag =_comon_tag(\@comonTagList, $feature);

				if($comonTag){ #Previous Level up feature had a comon tag
			
		    		push( @{$infoSequential->{$comonTag}{$id}{'level2'}}, $feature );			
			    	return $comonTag, $last_l1_f, $feature;								#### STOP HERE AND RETURN
				}
				else{
					######################
					# NEED THE LEVEL1 ID #
					my $l1_id="";			
					if($last_l1_f){ # case where previous level1 exists
						$l1_id=$last_l1_f->_tag_value('ID');
					}
					else{ # case where No level2 feature defined yet - I will need a bucketL2
						$l1_id = _create_ID($miscCount, $uniqID, $feature, $id, "nbis_noL1id");					
						$last_l1_f = clone($feature);
						create_or_replace_tag($last_l1_f,'ID',$l1_id); #modify Parent To keep only one
						$last_l1_f->primary_tag('gene');
					}
					my $l1_ID = lc($last_l1_f->_tag_value('ID'));
					push( @{$infoSequential->{$l1_ID}{$id}{'level2'}}, $feature );
					return $l1_ID , $last_l1_f, $feature;								#### STOP HERE AND RETURN
				}
			}

			#####################
	    	# Ckeck duplication #
    		if(! _it_is_duplication($duplicate, $omniscient, $uniqID, $feature)){

				############################################
				# keep track of link between level2->leve1 #
	  			if (! exists ($mRNAGeneLink->{lc($id)})){ 
					$mRNAGeneLink->{lc($id)}=lc($parent);
		 		}
		 		
		 		####################
		 		# SAVE THE FEATURE #
	      		push (@{$omniscient->{"level2"}{$primary_tag}{lc($parent)}}, $feature);
	      	}
      	}

      	###########
      	# NONE OF LEVEL FEATURE DEFINED
      	else{
      		print "gff3 reader warning: $primary_tag still not taken in account ! Please modify the code to define on of the three level of this feature.\n";
      	}	

    return $comonTag, $last_l1_f, $last_l2_f;
}

##==============================================================

sub _comon_tag{
	my ($comonTagList, $feature)=@_;

	my $result=undef;
   	foreach my $tag (@$comonTagList){

		if($feature->has_tag($tag)){
		    $result=lc($feature->_tag_value($tag));
		    return $result;
		}
	}
	return $result;
}

#feature is not yet saved in omniscient !
sub _it_is_duplication{
	my ($duplicate, $omniscient, $uniqID, $feature)=@_;

	my $is_dupli=undef;
	my $potentialList=undef;

	my $level = _get_level($feature);
	my $primary_tag = lc($feature->primary_tag);
	my $id = $uniqID->{$feature->_tag_value('ID')}; # check the original ID

	if(! exists_keys($omniscient,($level, $primary_tag, $id))){
	     return $is_dupli;
	}
	else{
		$potentialList=$omniscient->{$level}{$primary_tag}{$id}; #could be a feature if level1 or a list if level2/level3
	}

	#check if list is a feature or a reference of a list of feature
	my @list_feature; # will be a list of feature.
	if(ref($potentialList) eq 'ARRAY'){
		@list_feature=@$potentialList; #it's a reference of a list of feature
	}
	else{push (@list_feature, $potentialList);} # it's a feature

	#### PREPARE THE SENTENCE TO CHECK
	my $current_string=_create_string($uniqID, $feature);

	#Check all the level2 list element
	foreach my $feature_in_omniscient ( @list_feature ){

		my $string=_create_string($feature_in_omniscient);
		if($current_string eq $string){
			$is_dupli=1;
			push (@{$duplicate->{$level}{$primary_tag}{$id}}, $feature);
			delete $uniqID->{$feature->_tag_value('ID')}; # delete uniq ID that has been created for nothing
			last;
		}
	}

	if(! $is_dupli and $level eq "level1" and $omniscient->{"level1"}{$primary_tag}{$id}){
		warn "WARNING gff3 reader level1 : This feature level1 is not a duplicate (different than others feature level1) but has an ID already used. We cannot deal with that. @ the feature is: ".$feature->gff_string()."\n";  #Indeed If we change the ID we do not know wich sub-feature parent attribute value to modify (It could occur that we link sub feature not related)
	}
	return $is_dupli;
}

sub _get_level{
	my ($feature)=@_;

	my $primary_tag = lc($feature->primary_tag);

	if (exists(LEVEL1->{$primary_tag}) ) {
		return 'level1';
	}
	elsif(exists(LEVEL2->{$primary_tag}) ){
		return 'level2';
	}
	elsif(exists(LEVEL3->{$primary_tag}) ){
		return 'level3';
	}
	else{
		print "Error - cannot get the level of this feature !!";exit;
	}
}

sub _create_string{
	my ($uniqID, $feature)=@_;

	my $string=$feature->seq_id().$feature->start().$feature->end();

	my $primary_tag = lc($feature->primary_tag);
	if ( exists(LEVEL1->{$primary_tag}) ){
		$string .= $uniqID->{$feature->_tag_value('ID')}; # compare with original ID
	}
	if ( exists(LEVEL2->{$primary_tag}) ){
		$string .= $uniqID->{$feature->_tag_value('ID')}; # compare with original ID
		$string .= $uniqID->{$feature->_tag_value('Parent')}; # compare with original Parent
	}
	if ( exists(LEVEL3->{$primary_tag}) ){
		$string .= $uniqID->{$feature->_tag_value('Parent')}; # compare with original Parent
	}

	return $string;
}

# create an ID uniq. Don't give muilti-parent feature !
# If we have to create new ID for SPREADFEATURES they will not have a shared ID.
sub _check_uniq_id{
	my	($miscCount, $uniqID, $feature)=@_;
	
	my $id=undef;
	my $uID=undef;
	my $primary_tag = lc($feature->primary_tag);

	if($feature->has_tag('ID')){ #has the tag
		$id = $feature->_tag_value('ID');
		if(! SPREADFEATURE->{$primary_tag}){ #avoid CDS and UTR that can share identical IDs
			$uID = _create_ID($miscCount, $uniqID, $feature, $id, 'nbis_NEW')
		}
		else{$uID = $id;}
	}
	else{ #tag absent
		my $level = _get_level($feature);
		warn "gff3 reader error ".$level .": No ID attribute found @ for the feature: ".$feature->gff_string()."\n";
		$miscCount->{$primary_tag}++;
		$id = $primary_tag."-".$miscCount->{$primary_tag}; # create an ID and then check if not already taken
		$uID = _create_ID($miscCount, $uniqID, $feature, $id, 'nbis_NEW')
	}

	if(	$id ne $uID ){ #push the new ID if there is one
		create_or_replace_tag($feature, 'ID', $uID);
	}	
	#push the new ID	
	$uniqID->{$uID}=$id;

	return $uID;
}

# create the ID and add it to the feature.
sub _create_ID{
	my	($miscCount, $uniqID, $feature, $id, $prefix)=@_;
	
	my $key;

	if($prefix){
		$key=$prefix."-".lc($feature->primary_tag);
	}
	else{
		$key=lc($feature->primary_tag);
	}

	my $uID=$id;
	while(exists_keys($uniqID, $uID)){	 #loop until we found an uniq tag	
		$miscCount->{$key}++;
		$uID = $key."-".$miscCount->{$key};
	}
	return $uID;
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
	my ($hash_omniscient, $miscCount, $uniqID)=@_;

 	foreach my $tag_l1 (keys %{$hash_omniscient->{'level1'}}){
 	  	foreach my $id_l1 (keys %{$hash_omniscient->{'level1'}{$tag_l1}}){	     	
 		    my $neverfound="yes";
 		    foreach my $tag_l2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
 		        if ( exists ($hash_omniscient->{'level2'}{$tag_l2}{$id_l1} ) ){
 		          $neverfound=undef;last
 		        }   
 		    }
 		    if($neverfound){
 		    	warn "WARNING gff3 reader level1 : No child feature found for the $tag_l1 with ID @ ".$id_l1.".\n";
		    	
 		    	#check refseq case // Follw the gff3 standard but gap for feature L2 
     			check_refseq($hash_omniscient, $hash_omniscient->{'level1'}{$tag_l1}{$id_l1}, $miscCount, $uniqID);
 			    delete $hash_omniscient->{'level1'}{$tag_l1}{$id_l1}; # delete level1 // In case of refseq the thin has been cloned and modified, it is why we nevertheless remove it
		    }
 	 	}
 	}
}

# # All Level1 feature are for sure in omniscient, we have only the ID in sequential, other level feature are in sequential only if no parent has been found
# sub check_sequential{ # parcour de L3 to l1
# 	my ($infoSequential, $omniscient, $miscCount) = @_;
# 	print "check_sequential\n";
# 	foreach my $comonTag (keys %{$infoSequential} ){

# 		foreach my $bucket (keys %{$infoSequential->{$comonTag} } ){
# 			print "middle $comonTag $bucket\n";
# 			if ($bucket eq 'level1'){next;} #skip case level1 - strucutre of the hash different

# 			my $must_create_l2=undef;
# 			my $feature_L2 = undef;
# 			my $l2_ID=undef;
			
# 			if(! exists_keys($infoSequential,($comonTag, $bucket,'level3'))){
# 				print "Not normal, we had feature L1 or L2  without L3 feature associated. We skip it.\n";
# 				next;
# 			}
# 			else{
# 				foreach my $feature_L3 (@{$infoSequential->{$comonTag}{$bucket}{'level3'}} ){

# 					if(! exists_keys($infoSequential,($comonTag,$bucket,'level2'))  ){
# 							$must_create_l2=1;
# 							$feature_L2 = clone($infoSequential->{$comonTag}{$bucket}{'level3'}[0]);#create a copy of the first mRNA feature;
# 							$miscCount->{'noParent'}{'level2'}++;
# 							$l2_ID = 'level2-'.$miscCount->{'noParent'}{'level2'};
# 							create_or_replace_tag($feature_L2,'ID', $l2_ID); #modify ID to replace by parent value
# 					}
# 					else{
# 						$feature_L2 = $infoSequential->{$comonTag}{$bucket}{'level2'};
# 						$l2_ID = $feature_L2->_tag_value('ID');
# 					}
										
# 					my $primary_tag_L3 =  lc($feature_L3->primary_tag);
# 					create_or_replace_tag($feature_L3,'Parent', $l2_ID); #modify ID to replace by parent value
# 					push (@{$omniscient->{"level3"}{$primary_tag_L3}{lc($l2_ID)}}, $feature_L3);
# 				}
# 			}

# 			##
# 			# TAKE CARE OF LEVEL2
# 			#we necesserely have the level2 feature in sequential
# 			my $parentID = undef;
# 			if(! exists_keys($infoSequential,($comonTag,'level1'))){ #In that case level1 feature doesn't exists in $infoSequential and in $omniscient. I will be created by teh method check_gene_link_to_mrna
# 				$miscCount->{'noParent'}{'level1'}++;
# 				$parentID = "level1-".$miscCount->{'noParent'}{'level1'};
# 				print "pioioi\n";
# 			}
# 			else{	
# 				$parentID = $infoSequential->{$comonTag}{'level1'};
# 				print "My parent = $parentID \n";
# 			}

# 			if($must_create_l2){
# 				create_or_replace_tag($feature_L2,'Parent', $parentID ); # change parentID 
# 				if (! exists_keys($omniscient,("level3",'cds', lc($l2_ID)) )  ){
# 					$feature_L2->primary_tag('mRNA'); # change primary tag
# 				}
# 				else{	
# 					$feature_L2->primary_tag('RNA'); # NO CDS we cannot guess what kind of level2 feature it is. So we use this general term to describe it
# 				}
# 			}
# 			else{ #It exists, we use it as it is
# 				$feature_L2 = $infoSequential->{$comonTag}{$bucket}{'level2'};
# 			}

# 			my $primary_tag_L2 =  lc($feature_L2->primary_tag);	

# 			check_level2_positions($omniscient, $feature_L2);
# 			push (@{$omniscient->{"level2"}{$primary_tag_L2}{lc($parentID)}}, $feature_L2);
# 			print "end check_sequential\n";
# 		}
# 		#LEVEL 1 IS taking care later
# 	}

# }

# Level3 related to level1 wihtout level2 defined
# so the id of L1 is transferred to l2 to keep the parent attribute of l3 correct
sub check_refseq{ # refseq case as example
 	my ($omniscient, $l1_feature, $miscCount, $uniqID)=@_;

 	my $l1_ID=$l1_feature->_tag_value('ID');

 	foreach my $tag_l3 (keys %{$omniscient->{'level3'}}){
 	  	foreach my $id_l2 (keys %{$omniscient->{'level3'}{$tag_l3}}){
 	        if ($id_l2 eq lc($l1_ID) ){
		          
 		          	my $l2_feature = clone($l1_feature);#create a copy of the first mRNA feature;
 		          	if (exists_keys($omniscient,("level3",'cds', lc($l1_ID)) )  ){
 	               		$l2_feature->primary_tag('mRNA');
 	               	}
 	               	else{ #we cannot guess
 	               		$l2_feature->primary_tag('RNA');
 	               	}

 			        #modify Level1:
 			        $l1_ID = _create_ID($miscCount, $uniqID, $l1_feature, $l1_ID, 'nbis_NEW');
 				  	create_or_replace_tag($l1_feature,'ID', $l1_ID); #modify ID to replace by parent value
 				  	my $l1_feature_clone = clone($l1_feature);#create a copy of the first mRNA feature;
 				  	$omniscient->{"level1"}{lc($l1_feature->primary_tag)}{lc($l1_ID)} = $l1_feature_clone;

 				  	#Modify parent L2
 				  	create_or_replace_tag($l2_feature,'Parent', $l1_ID); #modify ID to replace by parent value
 				  	push(@{$omniscient->{"level2"}{lc($l2_feature->primary_tag)}{lc($l1_ID)}}, $l2_feature);
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
1;
