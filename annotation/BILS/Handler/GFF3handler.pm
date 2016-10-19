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

	A library to convert handle gff3 file and save it in memory 
	Inherits from 

=cut	

##########################
#### DEFINE CONSTANT #####
use constant LEVEL1 => { "gene" => 1, "sts" => 2, "match" => 3, "pseudogene" => 4,  "processed_pseudogene" => 4 };
use constant LEVEL2 => { "mrna" => 1, "ncrna" => 2, "mirna" => 3, "lcrna" => 4, "rrna" => 5, "srp_rna" => 6, "snrna" => 7, "lincrna" => 8, "trna" => 9, "trna_pseudogene" => 10, "snorna" => 11, "misc_rna" => 12, "rnase_p_rna" => 13, "tmrna" => 14, "match_part" => 15, "similarity" => 16, "rna" => 17, "pseudogenic_transcript" => 18, "transcript" => 19, "processed_transcript" => 20, "nmd_transcript_variant" => 21};
use constant LEVEL3 => { "cds" => 1, "exon" => 2, "stop_codon" => 3, "start_codon" => 4, "three_prime_utr" => 5, "five_prime_utr" => 6, "utr" => 7, "selenocysteine" => 8, "non_canonical_three_prime_splice_site" => 8, "non_canonical_five_prime_splice_site" => 10,
						"stop_codon_read_through" => 11, "sig_peptide" => 12, "tss" => 13, "tts" => 14, "intron" => 15 };
use constant SPREADFEATURE => {"cds" => 1, "three_prime_utr" => 2, "five_prime_utr" => 3, "utr" => 4};

# Save in omniscient hash (sorted in a specific way (3 levels)) a whole gff3 file
# Parser phylosophy: Parse by Parent/child ELSE 
#						Parse by comon_tag  ELSE
#							Parse by sequential (mean group feattures in a bucket, and the bucket change at each level2 feature, and bucket are join in a comon tag at each new L1 feature)
#Priority Parent > locus_tag.
sub slurp_gff3_file_JD {
	
	my ($self, $file, $comonTagAttribute, $verbose) = @_  ;
	
	my $gffio = Bio::Tools::GFF->new(-file => $file, -gff_version => 3);	


	### Handle to not print to much warning
	my %WARNS;
	my $nbWarnLimit=5;
  	local $SIG{__WARN__} = sub {
    my $message = shift;
    my @thematic=split /@/,$message ;
    
    $WARNS{$thematic[0]}++;
	    if($verbose){
	    	if ($WARNS{$thematic[0]} <= $nbWarnLimit){
				print $message;
			}
			if($WARNS{$thematic[0]} == $nbWarnLimit){
				print "$thematic[0] ************** Too much WARNING message we skip the next **************\n";
			}
	  	}
  	};
  	########### END WARNING LOCAL METHOD


	my %mRNAGeneLink; #Hast that keep track about link between l2 and l1
	my %omniscient; #Hast where all the feature will be saved
	my %duplicate;# Hash to store any counter. Will be use to create a new uniq ID 
	my %miscCount;# Hash to store any counter. Will be use to create a new uniq ID
	my %uniqID;# Hash to follow up with an uniq identifier every feature
	my %locusTAG;
	my %infoSequential;# Hash to store any counter. Will be use to create a new uniq ID
	my $locusTAGvalue=undef;
	my $last_l1_f=undef;
	my $last_l2_f=undef;
	my $last_l3_f=undef;

	#read every lines
	while( my $feature = $gffio->next_feature()) {
		($locusTAGvalue, $last_l1_f, $last_l2_f, $last_l3_f) = manage_one_feature($feature, \%omniscient, \%mRNAGeneLink, \%duplicate, \%miscCount, \%uniqID, \%locusTAG, \%infoSequential, $comonTagAttribute, $locusTAGvalue, $last_l1_f, $last_l2_f, $last_l3_f, 1);
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
    print "parsing finished - Start check\n";

    #Check sequential if we can fix cases. Hash to be done first, else is risky that we remove orphan L1 feature ... that are not yet linked to a sequential bucket
	if( keys %infoSequential ){ #hash is not empty
    	check_sequential(\%infoSequential, \%omniscient, \%miscCount, \%uniqID, \%mRNAGeneLink);
    }
    print "next1\n\n";
    #check level1 has subfeature else we remove it
  	remove_orphan_l1(\%omniscient, \%miscCount, \%uniqID, \%mRNAGeneLink); #or fix if level2 is missing (refseq case)
print "next2\n\n";
    #Check relationship between l3 and l2
    check_l3_link_to_l2(\%omniscient, \%mRNAGeneLink, \%miscCount, \%uniqID); # When creating L2 missing we create as well L1 if missing too
print "next3\n\n";
    #Check relationship L3 feature, exons have to be defined...
    check_exons(\%omniscient, \%miscCount, \%uniqID);
print "next4\n\n";
    #Check relationship between mRNA and gene.
    check_gene_link_to_mrna(\%omniscient);

print "next5\n\n";

    # To keep track of How many Warnings we got....
    foreach my $thematic (keys %WARNS){
  		my $nbW = $WARNS{$thematic};
  		if($nbW > $nbWarnLimit){
  			print "$nbW warning messages: $thematic\n";
  		}	
  	}

    #return
	return \%omniscient, \%mRNAGeneLink	;
}

#sub create_omniscient_from_feature_list {
#	my ($ref_featureList)=@_;
#
#	my %mRNAGeneLink;
#	my %omniscient;
#
#	foreach my $feature (@{$ref_featureList}) {
#		manage_one_feature($feature, \%omniscient, \%mRNAGeneLink);
 #   }
#
#	return \%omniscient, \%mRNAGeneLink	;
#}

# The method read a gff3 feature (one line), Check for correctly formated gff3 feature (ID and Parent attribute) 
sub manage_one_feature{
	
	my ($feature, $omniscient, $mRNAGeneLink, $duplicate, $miscCount, $uniqID, $locusTAG_uniq, $infoSequential, $comonTagAttribute, $last_locusTAGvalue, $last_l1_f, $last_l2_f, $last_l3_f, $debug)=@_;
		#print "Info before start: last_locusTAGvalue = $last_locusTAGvalue -- last_l2_f = $last_l2_f\n";
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
		if ($comonTagAttribute){push @comonTagList, $comonTagAttribute; } #add a new comon tag to the list if provided.
		my $locusTAGvalue=undef;

		####################################################
	    ########## Manage feature WITHOUT parent ###########	== LEVEL1 ==
	    ####################################################
	    if (exists(LEVEL1->{$primary_tag}) ) {
	    	if($debug){print "\nLEVEL1 case\n";}
	    	################
	    	# REPEATS case #
	    	if($primary_tag eq "match"){ #Manage REPEATS // Protein2genome // stuff in match/matchpart ###########	/!\ NO level3 features /!\ Compare to gene stuff we replace primary_tag by source_tag when filling the omniscient hash
	    		$primary_tag=$source_tag;
	    	}

	    	##########
			# get ID #
			#my $had_ID=undef;
	    	#if($feature->has_tag('ID')){$had_ID;}
	    	$id = lc(_check_uniq_id($miscCount, $uniqID, $feature));

	    	#####################
	    	# Ckeck duplication #
    		if(! _it_is_duplication($duplicate, $omniscient, $uniqID, $feature)){
    			
	    		################
	    		# Save feature #
	    		$last_l1_f = $feature;
	    		if($debug){print "Push-L1-omniscient level1 || $primary_tag || $id = feature\n";}
    			$omniscient->{"level1"}{$primary_tag}{$id}=$feature;
 	       		$locusTAG_uniq->{$id}=$id;

 	       		#############
 	       		# COMON TAG #
		        $locusTAGvalue =_get_comon_tag_value(\@comonTagList, $feature, $locusTAG_uniq);

				if($locusTAGvalue){
					if($debug){print "Push-L1-sequential $locusTAGvalue || level1 == $id\n";}
					$locusTAG_uniq->{$locusTAGvalue}=$id;
		    		$infoSequential->{$id}{'level1'}=$id;
		    	}
		    	else{
		    		$locusTAGvalue=$id;
		    	}

		    	#################
		    	#reinitialization
		    	$last_l2_f=undef; #opposite to I have a comon tag
		    	$last_l3_f=undef;
		    }

		    return $id, $last_l1_f, $last_l2_f, $last_l3_f;
	    }

      	###################################################
      	########## Manage feature WITHout CHILD  ##########		== LEVEL3 ==
      	###################################################
      	elsif ( exists(LEVEL3->{$primary_tag}) ){
      		if($debug){print "\nLEVEL3 case\n";}
      		##########
			# get ID #
	    	$id = _check_uniq_id($miscCount, $uniqID, $feature);
	    	print "id= $id || ".$feature->gff_string."\n";
	    	##############
			# get Parent #
			my @parentList;
      		if($feature->has_tag('Parent')){
      			@parentList = $feature->get_tag_values('Parent');
      			$locusTAGvalue=$last_locusTAGvalue;
			}
			else{ # In that case we create a uniq parentID to create a proper omniscient structure. But the feature itself stay intact without parentID.
				warn "WARNING gff3 reader level3: No Parent attribute found @ for the feature: ".$feature->gff_string()."\n";
				
				#################
				# COMON TAG PART1 
				$locusTAGvalue =_get_comon_tag_value(\@comonTagList, $feature, $locusTAG_uniq);

				######################
				# NEED THE LEVEL2 ID #
				my $l2_id="";
				# case where No level2 feature defined yet - I will need a bucketL2 
				# OR comon tag changed (= level1/level2 different) so we have to create a new level2 tag - but only if the last_comon tag is different as the parent of the last_l2_f (In that case we can use the last L2 feature. It was missing the comon tag in it). 
				if(! $last_l2_f or ($locusTAGvalue and ($locusTAGvalue ne $last_locusTAGvalue) and $last_locusTAGvalue ne lc($last_l2_f->_tag_value('Parent')) ) ){
					if($debug){print "Create L2 feature $locusTAGvalue ne $last_locusTAGvalue !\n";}
					# In case where no parent, no comon tag, and no sequential, we cannot deal at all with it !!!!
					if(! _get_comon_tag_value(\@comonTagList, $feature, $locusTAG_uniq) and ! $last_l1_f ){
						 print "Are you kiding ? ".$feature->gff_string()." is not a gff3 feature => BIG mistake ! this file cannot be parsed correclty. Eihter sequentialy or using comon tag.\n You could try to provide a correct comon tag (default is 'locus_tag')\n";
					}

					$l2_id = _create_ID($miscCount, $uniqID, $primary_tag, $id, "nbis_noL2id");					
					$last_l2_f = clone($feature);
					create_or_replace_tag($last_l2_f,'ID',$l2_id); #modify Parent To keep only one
					$last_l2_f->primary_tag('RNA');
				}
				else{ # case where previous level2 exists
					$l2_id=$last_l2_f->_tag_value('ID');
				}
				create_or_replace_tag($feature,'Parent',$l2_id); #modify Parent To keep only one

				#############
 	       		# COMON TAG  Part2
				if($locusTAGvalue){ #Previous Level up feature had a comon tag					
					if($debug){print "Push-L3-sequential-1 $locusTAGvalue || ".lc($l2_id)." || level3 == ".$feature->gff_string."\n";}
		    		push( @{$infoSequential->{$locusTAGvalue}{lc($l2_id)}{'level3'}}, $feature );
			    	return $locusTAGvalue, $last_l1_f, $last_l2_f, $feature;								#### STOP HERE AND RETURN
				}
				else{# No comon tag found 
					######################
					# NEED THE LEVEL1 ID #
					if(!$last_l1_f and $last_l3_f){ #particular case : Two l3 that follow each other, but first one has locus_tag but not the second
						if($debug){print "Push-L3-sequential-2 $last_locusTAGvalue || ".lc($l2_id)." || level3 == ".$feature->gff_string."\n";}
						push( @{$infoSequential->{$last_locusTAGvalue}{lc($l2_id)}{'level3'}}, $feature );
						return $last_locusTAGvalue, $last_l1_f, $last_l2_f, $feature;			
					}
					else{
						my $l1_id="";			
						if($last_l1_f){ # case where previous level1 exists
							$l1_id=$last_l1_f->_tag_value('ID');
						}
						else{ # case where No level1 feature defined yet - I will need a bucketL1
							$l1_id = _create_ID($miscCount, $uniqID, $primary_tag, $id, "nbis_noL1id");					
							$last_l1_f = clone($feature);
							create_or_replace_tag($last_l1_f,'ID',$l1_id); #modify Parent To keep only one
							$last_l1_f->primary_tag('gene');
						}

						#if($debug){print "Push-L3-Sequential-3 ".lc($l1_id)." || ".lc($l2_id)." || level3 == ".$feature->gff_string."\n";}
						#push( @{$infoSequential->{lc($l1_id)}{lc($l2_id)}{'level3'}}, $feature );
						if($debug){print "Push-L3-omiscient-3: level3 ".$primary_tag." || ".lc($l2_id)." == ".$feature->gff_string."\n";}
						push (@{$omniscient->{"level3"}{$primary_tag}{lc($l2_id)}}, $feature);
						return $l2_id, $last_l1_f, $last_l2_f, $feature;								#### STOP HERE AND RETURN
					}
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
							# if( ($primary_tag ne "cds") and (index($primary_tag, 'utr') == -1) ){ #avoid case where we can have the same ID (multi-feature)
							# 	my $substr = $feature->_tag_value('Parent');
							# 	my $clone_id=$feature->_tag_value('ID');
							# 	$clone_id=~ s/$substr//;
							# 	if($clone_id eq $feature->_tag_value('ID')){#Substring didn't work
							# 		$clone_id=$feature->_tag_value('ID')."$cptParent";
							# 	}else{$clone_id="$parent$clone_id";}
							# 	create_or_replace_tag($feature_clone,'ID',$clone_id); #modify Parent To keep only one
							# }

							_check_uniq_id($miscCount, $uniqID, $feature_clone); #Will change the ID if needed

							if($debug){print "Push-L3-omniscient-4 level3 || $primary_tag || ".lc($parent)." == feature_clone\n";}
							push (@{$omniscient->{"level3"}{$primary_tag}{lc($parent)}}, $feature_clone);
						}
						
						# It is the first parent. Do not clone the feature
						else{
							create_or_replace_tag($feature,'Parent',$parent); #modify Parent To keep only one
							if($debug){print "Push-L3-omniscient-5 level3 || $primary_tag || ".lc($parent)." == feature\n";}
							push (@{$omniscient->{"level3"}{$primary_tag}{lc($parent)}}, $feature);
						}
					}
					else{ #the simpliest case. One parent only
						if($debug){print "Push-L3-omniscient-6 level3 || $primary_tag || ".lc($parent)." == feature\n";}
						push (@{$omniscient->{"level3"}{$primary_tag}{lc($parent)}}, $feature);
					}
				}
				#Level3 key exists
				else{  # If not the first feature level3 with this primary_tag linked to the level2 feature
					# check among list of feature level3 already exits with an identical ID.
					# my $is_dupli=undef;
					# foreach my $feat ( @{ $omniscient->{"level3"}{$primary_tag}{lc($parent)} } ) {
					# 	# case where same feature is spread on different location (i.e utr, cds). In that case, following the gff3 gene ontologie specification, the ID can be share by all "peace of feature" building the feature entity.
					# 	if( _it_is_duplication($duplicate, $omniscient, $uniqID, $feature) ){
					# 			print "it is duplicated\n";
					# 			$is_dupli=1;
					# 			last;
					# 	}# case where a feature are an uniq element (exon, etc.)
					# 	elsif( lc($feat->_tag_value('ID')) eq lc($feature->_tag_value('ID')) ){
					# 		my $id = $feature->_tag_value('ID');

					# 		if( $feat->start() == $feature->start() and $feat->end() == $feature->end() ){
					# 			if($debug){print "Push-L3-duplicate-7 level3 || $primary_tag || ".lc($id)." == feature\n";}
					# 			push ( @{$duplicate->{"level3"}{$primary_tag}{lc($id)}}, $feature );
					# 			$is_dupli=1;
					# 			last;
					# 		}
					# 		else{ #case we have to change name ! GFF3 NON CONFORME
					# 			$miscCount->{$id}++;
					# 			$id=$id."-".$miscCount->{$id};
					# 			create_or_replace_tag($feature,'ID',$id); #modify Parent To keep only one
					# 		}
					# 	}
					# }
					#It is not a duplicated feature => save it in omniscient
					#if(! $is_dupli){
					if( ! _it_is_duplication($duplicate, $omniscient, $uniqID, $feature) ){
						# It is a multiple parent case
						if($allParent > 1){
							# Not the first parent, we have to clone the feature !!
							if($cptParent > 1){ 
								my $feature_clone=clone($feature);
								create_or_replace_tag($feature_clone,'Parent',$parent); #modify Parent To keep only one
								
								# #Take care of ID to have uniq one (case of exon,etc...)
								# if( ($primary_tag ne "cds") and (index($primary_tag, 'utr') == -1) ){
								# 	my $substr = $feature->_tag_value('Parent');
								# 	my $clone_id=$feature->_tag_value('ID');
								# 	$clone_id=~ s/$substr//;
								# 	if($clone_id eq $feature->_tag_value('ID')){#Substring didn't work
								# 		$clone_id=$feature->_tag_value('ID')."$cptParent";
								# 	}else{$clone_id="$parent$clone_id";}
								# 	create_or_replace_tag($feature_clone,'ID',$clone_id); #modify Parent To keep only one
								# }
								_check_uniq_id($miscCount, $uniqID, $feature_clone); #Will change the ID if needed

								if($debug){print "Push-L3-omniscient-8 level3 || $primary_tag || ".lc($parent)." == feature_clone\n";}
								push (@{$omniscient->{"level3"}{$primary_tag}{lc($parent)}}, $feature_clone);
							}
							# It is the first parent. Do not clone the feature
							else{
								create_or_replace_tag($feature,'Parent',$parent); #modify Parent To keep only one
								if($debug){print "Push-L3-omniscient-9 level3 || $primary_tag || ".lc($parent)." == feature\n";}
								push (@{$omniscient->{"level3"}{$primary_tag}{lc($parent)}}, $feature);
							}
						}
						else{ #the simpliest case. One parent only
							if($debug){print "Push-L3-omniscient-10 level3 || $primary_tag || ".lc($parent)." == feature\n";}
							push (@{$omniscient->{"level3"}{$primary_tag}{lc($parent)}}, $feature);
						}
					}
				}
      		}
      		return $last_locusTAGvalue, $last_l1_f, $last_l2_f, $feature;
      	}

      	##############################################
      	########## Manage feature the rest  ##########		== LEVEL2 ==
      	##############################################
      	elsif ( exists(LEVEL2->{$primary_tag}) ) {
    		if ($primary_tag eq "match_part" or $primary_tag eq "similarity") { ########## Manage REPEATS // Protein2genome // stuff in match/matchpart ###########	/!\ NO level3 features /!\ Compare to gene stuff we replace primary_tag by source_tag when filling the omniscient hash       	
				$primary_tag=$source_tag;  
    		}
    		if($debug){print "\nLEVEL2 case\n";}
    		#reinitialization
    		$last_l3_f=undef;

    		##########
			# get ID #
	    	$id = lc(_check_uniq_id($miscCount, $uniqID, $feature));

			##############
			# get Parent #
			if($feature->has_tag('Parent')){
				$parent = lc($feature->_tag_value('Parent'));
				$locusTAGvalue=$parent;
			}
			else{warn "WARNING gff3 reader level2 : No Parent attribute found for @ the feature: ".$feature->gff_string()."\n";

				#################
 	       		# COMON TAG PART1
				$locusTAGvalue =_get_comon_tag_value(\@comonTagList, $feature, $locusTAG_uniq);


				######################
				# NEED THE LEVEL1 ID #
				my $l1_ID="";
				# If I don't have a last_l1_f I create one. The Id can be used as comonTag. The feature can also be used later if comon_tag was existing, but mising for one of the feature. If we have comon tag, we check we are changing form the previous one before to create a new level1 feature. It's to deal with potential level2 (like mRNA isoforms).
				if(! $last_l1_f or ($locusTAGvalue and ($locusTAGvalue ne $last_locusTAGvalue) ) ){
					if($debug){print "create L1 feature\n";}
					$l1_ID = _create_ID($miscCount, $uniqID, $primary_tag, $id, "nbis_noL1id");					
					$last_l1_f = clone($feature);
					create_or_replace_tag($last_l1_f,'ID',$l1_ID); #modify Parent To keep only one
					$last_l1_f->primary_tag('gene');
				}
				else{ # case where previous level1 exists
					if($debug){print "take last L1 feature\n";}
					$l1_ID=$last_l1_f->_tag_value('ID');
				}
				create_or_replace_tag($feature,'Parent',$l1_ID); #modify Parent To keep only one

				#################
 	       		# COMON TAG PART2
				if($locusTAGvalue){ #Previous Level up feature had a comon tag
					if($debug){print "Push-L2-Sequential-1 $locusTAGvalue || ".lc($id)." || level2 == ".$feature->gff_string."\n";}
		    		$infoSequential->{$locusTAGvalue}{lc($id)}{'level2'} = $feature ;
		    		# keep track of link between level2->leve1 # Why using it ?
	  				#if (! exists ($mRNAGeneLink->{lc($id)})){ 
					#	$mRNAGeneLink->{lc($id)}=lc($locusTAGvalue);
		 			#}			
			    	return $locusTAGvalue, $last_l1_f, $feature, $last_l3_f;								#### STOP HERE AND RETURN
				}
				else{					
					#if($debug){print "Push-L2-Sequential-2 ".lc($l1_ID)." || ".lc($id)." || level2 == ".$feature->gff_string."\n";}
					#$infoSequential->{lc($l1_ID)}{lc($id)}{'level2'} = $feature ;
					if($debug){print "Push-L2-omniscient-2: level2 || ".$primary_tag." || ".lc($l1_ID)." == ".$feature->gff_string."\n";}
					push (@{$omniscient->{"level2"}{$primary_tag}{lc($l1_ID)}}, $feature);
					# keep track of link between level2->leve1 #
	  				if (! exists ($mRNAGeneLink->{lc($id)})){ 
						$mRNAGeneLink->{lc($id)}=lc($l1_ID);
		 			}	
					return lc($l1_ID) , $last_l1_f, $feature, $last_l3_f;								#### STOP HERE AND RETURN
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
		 		if($debug){print "Push-L2-omniscient-3 level2 || $primary_tag || $parent == feature\n";}
	      		push (@{$omniscient->{"level2"}{$primary_tag}{lc($parent)}}, $feature);
	      	}
	      	return $last_locusTAGvalue, $last_l1_f, $feature, $last_l3_f;
      	}

      	###########
      	# NONE OF LEVEL FEATURE DEFINED
      	else{
      		print "gff3 reader warning: $primary_tag still not taken in account ! Please modify the code to define on of the three level of this feature.\n";
      		return $locusTAGvalue, $last_l1_f, $last_l2_f, $last_l3_f;
      	}	

    print "Qrriver la cpas normal !!\n";exit;
}

##==============================================================

#check if the comom tag is present among the attributes
# return lower case value
sub _get_comon_tag_value{
	my ($comonTagList, $feature, $locusTAG_uniq)=@_;

	my $result=undef;
   	foreach my $tag (@$comonTagList){
		if($feature->has_tag($tag)){
		    $result=lc($feature->_tag_value($tag));
		    if(exists ($locusTAG_uniq->{$result}) ){
		    	return $locusTAG_uniq->{$result};
		    }
		    else{
		    	$locusTAG_uniq->{$result} = $result;
		    	return $result;
			}
		}
	}
	return $result;
}

#feature is not yet saved in omniscient !
sub _it_is_duplication{
	my ($duplicate, $omniscient, $uniqID, $feature)=@_;
	print $feature->gff_string."\n";
	my $is_dupli=undef;
	my $potentialList=undef;

	my $level = _get_level($feature);
	my $primary_tag = lc($feature->primary_tag);

	my $id = $uniqID->{$feature->_tag_value('ID')}; # check the original ID

	if($level eq "level1"){
		if(! exists_keys($omniscient,($level, $primary_tag, lc($id) ))){
			return $is_dupli; #return is not a dupli
		}
		else{
			$potentialList=$omniscient->{$level}{$primary_tag}{$id}; #push the feature L1 in potentialList
		}
	}
	else{ #feature l2 or l3

		my @parent = $feature->get_tag_values('Parent');
		foreach my $one_parent_ID (@parent){
			my $one_parent_uID = lc($uniqID->{$one_parent_ID}); # check the original ID

			foreach my $primary_tag ( keys %{$omniscient->{$level}} ){
				if (exists_keys($omniscient,($level, $primary_tag, $one_parent_uID))){
					push @{$potentialList}, @{$omniscient->{$level}{$primary_tag}{$one_parent_uID}};
				}
			}	
		}
		if(! $potentialList){ #potential list empty
		    return $is_dupli; #return is not a dupli
		}
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

		my $string=_create_string($uniqID, $feature_in_omniscient);
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

# find the level of the feature tested
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

#create string that should be uniq by feature
sub _create_string{
	my ($uniqID, $feature)=@_;

	my $string=$feature->seq_id().$feature->primary_tag().$feature->start().$feature->end();

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

# create an ID uniq. Don't give multi-parent feature !
# If we have to create new ID for SPREADFEATURES they will not have a shared ID.
sub _check_uniq_id{
	my	($miscCount, $uniqID, $feature)=@_;
	
	my $id=undef;
	my $uID=undef;
	my $primary_tag = lc($feature->primary_tag);

	if($feature->has_tag('ID')){ #has the tag
		$id = $feature->_tag_value('ID');
		if(! SPREADFEATURE->{$primary_tag}){ #avoid CDS and UTR that can share identical IDs
			$uID = _create_ID($miscCount, $uniqID, $primary_tag, $id, 'nbis_NEW'); #method will push the uID
			if(	$id ne $uID ){ #push the new ID if there is one
				create_or_replace_tag($feature, 'ID', $uID);
			}
			else{
				#push the uID	
				$uniqID->{$uID}=$id;
			}
		}
		else{
			$uID = $id;		
			#push the uID	
			$uniqID->{$uID}=$id;
		}

	}
	else{ #tag absent
		my $level = _get_level($feature);
		warn "gff3 reader error ".$level .": No ID attribute found @ for the feature: ".$feature->gff_string()."\n";
		$miscCount->{$primary_tag}++;
		$id = $primary_tag."-".$miscCount->{$primary_tag}; # create an ID and then check if not already taken
		$uID = _create_ID($miscCount, $uniqID, $primary_tag, $id, 'nbis_NEW'); #method will push the uID
		create_or_replace_tag($feature, 'ID', $uID);
	}

	return $uID;
}

# create the ID and add it to the feature.
sub _create_ID{
	my	($miscCount, $uniqID, $primary_tag, $id, $prefix)=@_;
	
	my $key;

	if($prefix){
		$key=$prefix."-".$primary_tag;
	}
	else{
		$key=$primary_tag;
	}

	my $uID=$id;
	while(exists_keys($uniqID, ($uID) )){	 #loop until we found an uniq tag	
		$miscCount->{$key}++;
		$uID = $key."-".$miscCount->{$key};
	}

	#push the new ID	
	$uniqID->{$uID}=$id;

	return $uID;
}

# check if mRNA have is Parental gene existing. If not we create it.
sub check_gene_link_to_mrna{
	my ($hash_omniscient)=@_;

	foreach my $primary_tag_l2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
		foreach my $id_l1 (keys %{$hash_omniscient->{'level2'}{$primary_tag_l2}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
			my $l1_exist=undef;
			foreach my $primary_tag_l1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
				if(exists_keys ($hash_omniscient, ('level1', $primary_tag_l1, $id_l1))){
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

# @Purpose: Remove the level1 feature that havn't subfeature linked to it. Before to remove it check if L3 is linked to it. In that case it is a format error that we will fix !
# @input: 1 => hash(omniscient hash)
# @output: none
sub remove_orphan_l1{
	my ($hash_omniscient, $miscCount, $uniqID, $mRNAGeneLink)=@_;

 	foreach my $tag_l1 (keys %{$hash_omniscient->{'level1'}}){
 	  	foreach my $id_l1 (keys %{$hash_omniscient->{'level1'}{$tag_l1}}){	     	
 		    my $neverfound="yes";
 		    foreach my $tag_l2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
 		        if ( exists_keys ( $hash_omniscient,('level2',$tag_l2,$id_l1) ) ){
 		          $neverfound=undef;last
 		        }   
 		    }
 		    if($neverfound){
 		    	
		    	
 		    	#check refseq case // Follw the gff3 standard but gap for feature L2 
     			if ( ! refseq_case($hash_omniscient, $hash_omniscient->{'level1'}{$tag_l1}{$id_l1}, $miscCount, $uniqID, $mRNAGeneLink)){
     				warn "WARNING gff3 reader level1 : No child feature found for the $tag_l1 with ID @ ".$id_l1.".\n";
     			}
 			    delete $hash_omniscient->{'level1'}{$tag_l1}{$id_l1}; # delete level1 // In case of refseq the thin has been cloned and modified, it is why we nevertheless remove it
		    }
 	 	}
 	}
}

# Level3 related to level1 wihtout level2 defined
# so the id of L1 is transferred to l2 to keep the parent attribute of l3 correct
sub refseq_case{ # refseq case as example
 	my ($omniscient, $l1_feature, $miscCount, $uniqID ,$mRNAGeneLink)=@_;

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
 			        $l1_ID = _create_ID($miscCount, $uniqID, lc($l1_feature->primary_tag), $l1_ID, 'nbis_NEW');
 				  	create_or_replace_tag($l1_feature,'ID', $l1_ID); #modify ID to replace by parent value
 				  	my $l1_feature_clone = clone($l1_feature);#create a copy of the first mRNA feature;
 				  	$omniscient->{"level1"}{lc($l1_feature->primary_tag)}{lc($l1_ID)} = $l1_feature_clone;

 				  	#Modify parent L2
 				  	create_or_replace_tag($l2_feature,'Parent', $l1_ID); #modify ID to replace by parent value
 				  	push(@{$omniscient->{"level2"}{lc($l2_feature->primary_tag)}{lc($l1_ID)}}, $l2_feature);

 				  	#fill the $mRNAGeneLink hash
 				  	$mRNAGeneLink->{lc($id_l2)} = lc($l1_ID); # Always need to keep track about l2->l1, else the method check_l3_link_to_l2 will recreate a l1 thinking this relationship is not fill

 				  	return 1;
 		    }
 	 	}
 	}
 	return 0;
}

# @Purpose: Check relationship betwwen L3 and L2. If L2 es missing we create it. When creating L2 missing we create as well L1 if missing too.
# @input: 4 => hash(omniscient hash), hash(mRNAGeneLink hash), hash(miscCount hash), hash(uniqID hash)
# @output: none
sub check_l3_link_to_l2{
	my ($hash_omniscient, $mRNAGeneLink, $miscCount, $uniqID)=@_;

 	foreach my $tag_l3 (keys %{$hash_omniscient->{'level3'}}){
 	  	foreach my $id_l2 (keys %{$hash_omniscient->{'level3'}{$tag_l3}}){

 	  		#check if L2 exits
 	  		if (! exists($mRNAGeneLink->{$id_l2}) ) {
 	  			
 	  		#start fill L2
 	  			my $l2_feature=clone($hash_omniscient->{'level3'}{$tag_l3}{$id_l2}[0]);#create a copy of the first mRNA feature;
				my $new_ID = $l2_feature->_tag_value('Parent');
				create_or_replace_tag($l2_feature,'ID', $new_ID); #modify ID to replace by parent value
				my $primary_tag_l2;
				if( exists_keys($hash_omniscient,('level3', 'cds', $id_l2)) ) {
					$primary_tag_l2="mRNA";
				}
				else{
					$primary_tag_l2="RNA";
				}
				$l2_feature->primary_tag($primary_tag_l2); # change primary tag
				check_level2_positions($hash_omniscient, $l2_feature);	# check start stop if isoforms exists

			#fill L1
				my $l1_feature=clone($hash_omniscient->{'level3'}{$tag_l3}{$id_l2}[0]);#create a copy of the first mRNA feature;
				$l1_feature->remove_tag('Parent'); # remove parent ID because, none.
				
				#Deal case where we reconstruct other thing than a gene
				my $primary_tag_l1=undef;
				if(lc($l1_feature->primary_tag) =~ /match/){ $primary_tag_l1="match"; }
				else{ $primary_tag_l1="gene"; }
				$l1_feature->primary_tag($primary_tag_l1); # change primary tag

				my $new_ID_l1 = _check_uniq_id($miscCount, $uniqID, $l1_feature);
				create_or_replace_tag($l1_feature,'ID', $new_ID_l1); #modify ID to replace by parent value
				
			#finish fill Level2
				create_or_replace_tag($l2_feature, 'Parent', $new_ID_l1); # remove parent ID because, none.
				#save new feature L2
				push (@{$hash_omniscient->{"level2"}{$primary_tag_l2}{lc($new_ID_l1)}}, $l2_feature);

			#finish fill Level1
				check_level1_positions($hash_omniscient, $l1_feature);	# check start stop if isoforms exists
				#save new feature L1
				$hash_omniscient->{"level1"}{$primary_tag_l1}{lc($new_ID_l1)} = $l1_feature; # now save it in omniscient

 	  		}
 	  	}
 	}	 
}

# @Purpose: Check L3 features. If exon are missing we create them.
# @input: 3 =>  hash(omniscient hash), hash(miscCount hash), hash(uniqID hash)
# @output: none
sub check_exons{
	my ($hash_omniscient, $miscCount, $uniqID)=@_;

	my %checked;
	foreach my $tag_l3 (keys %{$hash_omniscient->{'level3'}}){
		if ($tag_l3 ne "exon"){
 	  		foreach my $id_l2 (keys %{$hash_omniscient->{'level3'}{$tag_l3}}){
 	  			if( ! exists_keys(\%checked,($id_l2)) ){ #l2 already checked

 	  				#Case where No exon feature among L3
	 	  			#if( ! exists_keys($hash_omniscient,('level3','exon', $id_l2)) ){ #l3 feature but no exon among them... need to recreate them.

	 	  				my $feature_example=undef; # will be used to create the exon features
	 	  				my $list_location=[];

	 	  				foreach my $tag_l3 (keys %{$hash_omniscient->{'level3'}}){

				 	  		if( exists_keys($hash_omniscient,('level3',$tag_l3, $id_l2)) ){				 	  		

				 	  			foreach my $l3_feature (@{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2}}){

				 	  				if(! $feature_example){
				 	  					$feature_example=$l3_feature;
				 	  				}
				 	  				my $locationRefList=[[$l3_feature->start, $l3_feature->end]];
				 	  				$list_location = _manage_location($locationRefList, $list_location, 1);
				 	  			}
				 	  		}
				 	  	}
				 	  	$checked{$id_l2}++;
				 	  	print  "LEst continue: ".Dumper($list_location);exit;
				 	  	foreach my $location (@{$list_location}){
				 	  		my $feature_exon = clone($feature_example);#create a copy of a random feature l3;
							$feature_exon->start($location->[0]);
				 	  		$feature_exon->end($location->[1]);
				 	  		$feature_exon->primary_tag('exon');
				 	  		my $uID = _check_uniq_id($miscCount, $uniqID, $feature_exon);
				 	  		create_or_replace_tag($feature_exon, 'ID', $uID); # remove parent ID because, none.
							#save new feature L2
							push (@{$hash_omniscient->{"level3"}{'exon'}{$id_l2}}, $feature_exon);
				 	  	} 
	 	  		#	}
	 	  		}
 	  		}
 	  	}
 	}
}

# @Purpose: Will merge a list of "location" (tuple of integer), and another list of location. If two location overlap or are adjacent, pmly one location will be kept that represent the most extrem values
# @input: 3 =>  list of integer tuple ([[X,Y][Z,W]] or [[X,Y]]),  list of integer tuple, verbose option for debug
# @output: list of list of integer tuple
sub _manage_location{
	my ($locationRefList, $locationTargetList, $verbose) = @_;

	if($verbose){print "Enter Ref: ".Dumper($locationRefList)."\n"; print "Enter Target: ".Dumper($locationTargetList)."\n";}

	my @new_location_list; #new location list that will be returned once filled
	
	if (@$locationTargetList >= 1){ #check number of location -> List not empty
		
		my $check_list=1;

		while($check_list){
			$check_list=undef;
			
			my %tupleChecked=();
			foreach my $location_from_ref_list (@{$locationRefList}){

				my $location_skipped=undef; #keep track for later. If it has never been further, we will not save it later, because it will cause duplicates
				foreach my $location_from_target_list (@{$locationTargetList}){
					
					#skip particular cases						
						#case test againt himself
						if($location_from_ref_list->[0] ==  $location_from_target_list->[0] and $location_from_ref_list->[1] ==  $location_from_target_list->[1]){ # skip check again himself
							if($verbose){print "skip himself\n";}
							$location_skipped=1;
							next;
						}
						#case test already done
						my @tuple = sort {$a <=> $b} (@$location_from_ref_list,@$location_from_target_list);				
						my $tupleString = join(' ', @tuple);

						if(exists_keys (\%tupleChecked,($tupleString)) ) {
							if($verbose){print "skip tested\n";}
							$location_skipped=1;
							next;
						}
						else{
							$tupleChecked{$tupleString}++;
						} #track tuple checked when we conpare a lsit aginst itself. Because a test A<=>B is equal to B<=>A
	
					#check if we modify the location or not
					my ($new_location, $overlap) = _manage_location_lowLevel($location_from_ref_list, $location_from_target_list);
					if($verbose){print "push1\n";}
					push @new_location_list, [@$new_location];

					if( ($new_location->[0] !=  $location_from_target_list->[0] or $new_location->[1] !=  $location_from_target_list->[1]) or ($overlap)){ # location has been modified or not modifier but overlap (It means overlap completly ... it is included in)
						$check_list=1;
					}

					
				}
				if(! $check_list and ! $location_skipped){ #No modification done
					if($verbose){print "push0\n";}
					push @new_location_list, [@{$location_from_ref_list}];
				}
				#print "location_skipped $location_skipped\n";
			}
			if($check_list){ #modification done we have to re-check all location between each other
				
				if (scalar @new_location_list == 1){
					$check_list = undef; #Do not check the list if it is size one ! Because we will be stuck => case test againt himself is skipped !
				}
				else{
					$locationTargetList = [@new_location_list];
					$locationRefList = [@new_location_list];
					if($verbose){print "clear new_location_list\n";}
					@new_location_list=();
					%tupleChecked=();
				}
			}
		}
	}
	else{#check number of location -> none
		if($verbose){print "returnA: ".Dumper($locationRefList)."\n";}
		return $locationRefList;
	}
	if($verbose){print "returnB: ".Dumper(\@new_location_list)."\n";}
	return \@new_location_list;
}

# @Purpose: Modify the location2 if it overlap the location1 by keeping the extrem values. Return the location2 intact if no overlap.
# @input: 2 =>  integer tuple [X,Y],  list of integer tuple
# @output: 2 => ref of a list of 2 element, boolean
sub _manage_location_lowLevel{
	my ($location, $location2) = @_;

	my $new_location = [@{$location2}];
	my $overlap=undef;

	if ( ($location2->[0] <= $location->[1]+1) and ($location2->[1]+1 >= $location->[0]) ){ #it overlap or are consecutive

		$overlap=1;

		if($location2->[0] > $location->[0]){
			$new_location->[0]=$location->[0];
		}
		if($location->[1] > $location2->[1]){
			$new_location->[1]=$location->[1];
		}
	}

	return $new_location, $overlap;
}

# # All Level1 feature are for sure in omniscient, we have only the ID in sequential, other level feature are in sequential only if no parent has been found
sub check_sequential{ # Goes through from L3 to l1
 	my ($infoSequential, $omniscient, $miscCount, $uniqID, $mRNAGeneLink) = @_;
 	#print Dumper($infoSequential)."dumpermm\n";
 	#print Dumper($omniscient)."dumper omniscient\n";
 	foreach my $comonTag (keys %{$infoSequential} ){ #comon tag was l1 id wheb no real comon tag present
 		print "comonTag $comonTag\n";
 		foreach my $bucket (keys %{$infoSequential->{$comonTag} } ){ #bucket = level1 or Id L2
 			print "comonTag $comonTag bucket $bucket\n";
 			if ($bucket eq 'level1'){next;} #skip case level1 - structure of the hash different

 			my $must_create_l2=undef;
 			my $feature_l2 = undef;

			#Bucket is an uniq ID created during the reading process. So it can be used as uniq ID.
 			if(! exists_keys($infoSequential,($comonTag, $bucket, 'level3') ) ){
 				print "Not normal, we had feature L1 or L2  without L3 feature associated. We skip it.\n"; #We cannot guess the structure except if it is prokaryote... should we improve that ?
 				next;
 			}
			else{
 				foreach my $feature_L3 (@{$infoSequential->{$comonTag}{$bucket}{'level3'}} ){

 					if(! exists_keys($infoSequential,($comonTag, $bucket,'level2'))  ){
 						print "level2 does not exits in sequential !\n";

 						#take L2 from omniscient if already exits
 						if(exists($mRNAGeneLink->{lc($bucket)}) ){
 							print "rien \n";
 							my $l1_id = $mRNAGeneLink->{lc($bucket)};
 							foreach my $tag_l2 (keys %{$omniscient->{'level2'}}){
 								if(exists_keys($omniscient, ('level2', $tag_l2, lc($l1_id) ) ) ){
		 							foreach my $featureL2 (@{$omniscient->{'level2'}{$tag_l2}{lc($l1_id)}}){
		 								if(lc($featureL2->_tag_value('ID')) eq $bucket){
		 									print "level2 exits in omniscient !\n";
		 									$feature_l2 = $featureL2;
		 									last;
		 								}
		 							}
		 							if($feature_l2){last;}
		 						}
	 						}
 							$infoSequential->{$comonTag}{$bucket}{'level2'} = $feature_l2;
						}
 						else{#create l2
 							print "create level2  !\n";
	 						$must_create_l2=1;
	 						$feature_l2 = clone($infoSequential->{$comonTag}{$bucket}{'level3'}[0]);#create a copy of the first mRNA feature;
							#manage primary tag
							my $primary_tag_l2='RNA';
							foreach my $feature_L3 (@{$infoSequential->{$comonTag}{$bucket}{'level3'}} ){

	 							if ( lc($feature_L3->primary_tag) eq 'cds'){
	 								$primary_tag_l2 ='mRNA';
	 								last;
	 							}
	 						}
	 						$feature_l2->primary_tag($primary_tag_l2);

	 						#Manage ID 
								create_or_replace_tag($feature_l2,'ID', $bucket); #modify ID to replace by parent value
							#Manage Parent
								my $parentID = undef;
							 	if( exists_keys($infoSequential,($comonTag,'level1'))  ){
	 								$parentID = lc($infoSequential->{$comonTag}{'level1'});
	 							}
								else{
									my $IDgoodCast = _id_exists_in_l1_omniscient($omniscient, $comonTag);
									if($IDgoodCast){
											$parentID = $IDgoodCast;
									}

						 			if( ! $parentID ){ #In that case level1 feature doesn't exists in $infoSequential and in $omniscient. I will be created by teh method check_gene_link_to_mrna 
						 				#my	($miscCount, $uniqID, $primary_tag, $id, $prefix)=@_;
						 				$parentID =  _create_ID($miscCount, $uniqID, 'gene', "gene-1", 'nbis_NEW');
						 			}
						 		}
						 		print "Parent ID created = $parentID\n";
					 			create_or_replace_tag($feature_l2,'Parent', $parentID ); # change parentID 

					 		push (@{$omniscient->{"level2"}{lc($primary_tag_l2)}{lc($parentID)}}, $feature_l2);
					 		$mRNAGeneLink->{lc($bucket)} = lc($parentID); # Always need to keep track about l2->l1, else the method check_l3_link_to_l2 will recreate a l1 thinking this relationship is not fill
					 		$infoSequential->{$comonTag}{$bucket}{'level2'} = $feature_l2;
					 	}
 					}
					else{
						
						#MUST push L2 in omniscient if absent !
						$feature_l2=$infoSequential->{$comonTag}{$bucket}{'level2'};
						#print "level2 exits in sequential - $comonTag $bucket! ".$feature_l2->gff_string."\n";

						if(! exists($mRNAGeneLink->{lc($bucket)}) ){
							#print "level2 does not exits in mRNAGeneLink(omniscient) !".$feature_l2->gff_string."\n";
							push (@{$omniscient->{"level2"}{lc($feature_l2->primary_tag)}{lc($feature_l2->_tag_value('Parent'))} }, $feature_l2);
							$mRNAGeneLink->{lc($feature_l2->_tag_value('ID'))} = lc($feature_l2->_tag_value('Parent'));
						}
					}					
 					my $primary_tag_L3 =  lc($feature_L3->primary_tag);
 					create_or_replace_tag($feature_L3,'Parent', $feature_l2->_tag_value('ID')); #modify ID to replace by parent value

 					print "push-omniscient: level3 || ".$primary_tag_L3." || ".lc($bucket)." == ".$feature_L3->gff_string."\n";
 					push (@{$omniscient->{"level3"}{$primary_tag_L3}{lc($bucket)}}, $feature_L3);
 				}
 			}

 			if($must_create_l2){
 				check_level2_positions($omniscient, $feature_l2);
 			}
 		}
 		#LEVEL 1 IS taking care later
 	}
 	#print Dumper($omniscient)."dumper omniscient\n";
}

#print return ID if exists original cast
sub _id_exists_in_l1_omniscient{
	my ($omniscient, $id)=@_;

	my $id_good_cast=undef;
	foreach my $tag_l1 (keys %{$omniscient->{'level1'}} ){
		if(exists_keys($omniscient, ('level1',$tag_l1, $id))){
			$id_good_cast = $omniscient->{'level1'}{$tag_l1}{$id}->_tag_value('ID');
			return $id_good_cast;
		}
	}
	return $id_good_cast;
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
