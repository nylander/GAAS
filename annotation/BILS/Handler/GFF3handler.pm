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
our @EXPORT_OK   = qw(check_gene_positions check_mrna_positions modelate_utr_and_cds_features_from_exon_features_and_cds_start_stop create_omniscient slurp_gff3_file_JD create_omniscient_from_feature_list);
our %EXPORT_TAGS = ( DEFAULT => [qw()],
                 Ok    => [qw(check_gene_positions check_mrna_positions modelate_utr_and_cds_features_from_exon_features_and_cds_start_stop create_omniscient slurp_gff3_file_JD create_omniscient_from_feature_list)]);
=head1 SYNOPSIS



=head1 DESCRIPTION

	A library to convert 
	Inherits from 

	Dont take in account repeat and multi parent feature!!!
	
=cut	

# Save in omniscient hash (sorted in a specific way (3 levels)) a whole gff3 file
sub slurp_gff3_file_JD {
	
	my ($self, $file) = @_  ;
	
	my $gtfio = Bio::Tools::GFF->new(-file => $file, -gff_version => 3);	

	### Handle to not print to much warning
	my %WARNS;
	my $nbWarn=5;
  	local $SIG{__WARN__} = sub {
    my $message = shift;
    my @thematic=split /@/,$message ;
    
    $WARNS{$thematic[0]}++;
    	if ($WARNS{$thematic[0]} < $nbWarn){
			print $message;
		}
		elsif($WARNS{$thematic[0]} == $nbWarn){
			print "$thematic[0] ************** Too much WARNING message we skip the next **************\n";
		}
  	};

	my %mRNAGeneLink;
	my %omniscient;
	my %duplicate;# Hash to store any counter. Will be use to create a new uniq ID 
	my %miscCount;# Hash to store any counter. Will be use to create a new uniq ID



	while( my $feature = $gtfio->next_feature()) {
		manage_one_feature($feature, \%omniscient, \%mRNAGeneLink, \%duplicate, \%miscCount);
    }

    #Check if duplicate detected:
    my $keyExist = keys %duplicate;
    if($keyExist){#print result
    	my $stringPrint .= "################################################################\n# Achthung /\\ Attention /\\ Be carefull => ID duplicate found ! #\n".
      	"# Duplicated features have been removed (Keep only one per ID) #\n################################################################\n\n";
      	print $stringPrint;
      	my $gffout= Bio::Tools::GFF->new( -fh => \*STDOUT );
      	print_duplicates(\%duplicate, \%omniscient, $gffout);
      	print "We removed them from the regular output.\n";
    }

    #Check relationship between mRNA and gene
    check_gene_link_to_mrna(\%omniscient);


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
	
	my ($feature, $omniscient, $mRNAGeneLink, $duplicate, $miscCount)=@_;

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
	    	if($feature->has_tag('ID')){
	    		$id = lc($feature->_tag_value('ID')) ;
	    	}
	    	else{print "error !! No ID attribute found for the feature $feature->gff_string()\n";}
    	
    		#Save feature
    		if(! exists_keys($omniscient,('level1',$primary_tag,$id))){
	        	$omniscient->{"level1"}{$primary_tag}{$id}=$feature;
	        }
	        else{push (@{$duplicate->{"level1"}{$primary_tag}{$id}}, $feature);}
	    }

      	###################################################
      	########## Manage feature WITHout CHILD  ##########		== LEVEL3 ==
      	###################################################

      	elsif ( ($primary_tag eq "cds") or ($primary_tag eq "exon") or ($primary_tag eq "stop_codon") or ($primary_tag eq "start_codon") or ($primary_tag eq "three_prime_utr") or ($primary_tag eq "five_prime_utr") or ($primary_tag eq "utr"),
      	 or ($primary_tag eq "selenocysteine") or ($primary_tag eq "non_canonical_three_prime_splice_site") or ($primary_tag eq "non_canonical_five_prime_splice_site") or ($primary_tag eq "stop_codon_read_through"),
      	 or ($primary_tag eq "sig_peptide") ){

      		# manage ID
      		if(! $feature->has_tag('ID')){
      			warn "WARNING gff3 reader level3: No ID attribute found @ for the feature".$feature->gff_string()."\n";
				$miscCount->{'noID'}{$primary_tag}++;
				my $l3_id = $primary_tag."-".$miscCount->{'noID'}{$primary_tag};
				$feature->add_tag_value('ID', $l3_id);	
			}

			# manage Parent
			my @parentList;
      		if($feature->has_tag('Parent')){
      			@parentList = $feature->get_tag_values('Parent');
			}
			else{ # In that case we create a uniq parent ID to create a preper omniscient structure. But the feature itself stay intact without parentID.
				warn "WARNING gff3 reader level3: No Parent attribute found @ for the feature".$feature->gff_string()."\n";
				$miscCount->{'noParent'}{$primary_tag}++;
				$parent = $primary_tag."-".$miscCount->{'noParent'}{$primary_tag};
				push(@parentList, $parent);
			}

			#Save feature and check duplicates	(treat also cases where there is multiple parent. => In that case we expand to create a uniq feature for each)	
      		foreach my $parent (@parentList){ # first feature level3 with this primary_tag linked to the level2 feature
				if(! exists_keys($omniscient,('level3',$primary_tag,lc($parent)))){
					# save feature in omciscient
					create_or_replace_tag($feature,'Parent',$parent); #modify Parent To keep only one
					push (@{$omniscient->{"level3"}{$primary_tag}{lc($parent)}}, $feature);
				}
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
					if(! $is_dupli){
						# save feature in omniscient
						create_or_replace_tag($feature,'Parent',$parent); #modify Parent To keep only one
						push (@{$omniscient->{"level3"}{$primary_tag}{lc($parent)}}, $feature);
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

			else{warn "gff3 reader error level2: No ID attribute found @ for the feature ".$feature->gff_string()."\n";
				$miscCount->{'noID'}{$primary_tag}++;
				$id = $primary_tag."-".$miscCount->{'noID'}{$primary_tag};		
			}

			#get Parent
			if($feature->has_tag('Parent')){
				$parent = lc($feature->_tag_value('Parent'));
			}
			else{warn "WARNING gff3 reader level2 : No Parent attribute found for @ the feature".$feature->gff_string()."\n";
				$miscCount->{'noParent'}{$primary_tag}++;
				$parent = $primary_tag."-".$miscCount->{'noParent'}{$primary_tag};	
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
      			my $is_dupli=undef;
      			foreach my $feat ( @{$omniscient->{"level2"}{$primary_tag}{lc($parent)} }){
					if($feat->_tag_value('ID') eq $feature->_tag_value('ID')){
						push (@{$duplicate->{"level2"}{$primary_tag}{lc($parent)}}, $feature);
						my $is_dupli=1;
						last;
					}
     			}
     			if(! $is_dupli){ # No similar ID found, we keep it as isoform
	     			push (@{$omniscient->{"level2"}{$primary_tag}{lc($parent)}}, $feature);
	     		}		
      		}

      	}

      	#####################################
	    ########## Manage REPEATS ###########	/!\ NO level3 features /!\ Compare to gene stuff we replace primary_tag by source_tag when filling the omniscient hash
	    #####################################
      	elsif($source_tag =~ /repeat/) {

      		#									== LEVEL1 == 
			if ($primary_tag eq "match" ) {
	      		#get ID
		    	if($feature->has_tag('ID')){
		    		$id = lc($feature->_tag_value('ID')) ;
		    	}
		    	else{print "error !! No ID attribute found for the feature $feature->gff_string()\n";}
	    	
	    		#Save feature
	    		if(! exists_keys($omniscient,('level1',$source_tag,$id))){
		        	$omniscient->{"level1"}{$source_tag}{$id}=$feature;
		        }
		        else{push (@{$duplicate->{"level1"}{$source_tag}{$id}}, $feature);}
	    	}
	    	      		#									== LEVEL2 == 
			elsif ($primary_tag eq "match_part" ) {

				#get ID
		    	if($feature->has_tag('ID')){
			    	my @values = $feature->get_tag_values('ID');
					$id = lc(shift @values) ;
				}
				else{warn "gff3 reader error level2: No ID attribute found @ for the feature ".$feature->gff_string()."\n";
					$miscCount->{'noID'}{$source_tag}++;
					$id = $source_tag."-".$miscCount->{'noID'}{$source_tag};		
				}

				#get Parent
				if($feature->has_tag('Parent')){
					my @values = $feature->get_tag_values('Parent');
					$parent = lc(shift @values);
				}
				else{warn "WARNING gff3 reader level2 : No Parent attribute found for @ the feature".$feature->gff_string()."\n";
					$miscCount->{'noParent'}{$source_tag}++;
					$parent = $source_tag."-".$miscCount->{'noParent'}{$source_tag};	
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
	      			my $is_dupli=undef;
	      			foreach my $feat ( @{$omniscient->{"level2"}{$source_tag}{lc($parent)} }){
						if($feat->_tag_value('ID') eq $feature->_tag_value('ID')){ # on similar ID exits. We save it as a duplicate
							push (@{$duplicate->{"level2"}{$source_tag}{lc($parent)}}, $feature);
							my $is_dupli=1;
							last;
						}
	     			}
	     			if(! $is_dupli){ # No similar ID found, we keep it as isoform
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
				warn "WARNING gff3 reader level2 : No Parent feature found with the ID".$id_l1.". We will create one.\n";
				my $gene_feature=clone($hash_omniscient->{'level2'}{$primary_tag_l2}{$id_l1}[0]);#create a copy of the first mRNA feature;
				my $new_ID = $gene_feature->_tag_value('Parent');
				create_or_replace_tag($gene_feature,'ID', $new_ID); #modify ID to replace by parent value
				$gene_feature->remove_tag('Parent'); # remove parent ID because, none.
				check_level1_positions($hash_omniscient, $gene_feature);	# check start stop if isoforms exists
				$hash_omniscient->{"level1"}{'gene'}{lc($new_ID)}=$gene_feature; # now save it in omniscient
			}
		}
	}
}

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
 			$cds_feature->primary_tag('cds');
 			#get old name
 			my @IDs = $cds_feature->get_tag_values('ID');
 			my $ID = $IDs[0];
 			create_or_replace_tag($cds_feature,'ID',$ID.'-cds-'.$cds_counter); #modify name
 			push(@cds_features, $cds_feature);#save that cds
 			$cds_counter++;
 			if($exon_feature->start < $ORFstart){
 				my $utr_feature=clone($exon_feature);#create a copy of the feature 
 				$utr_feature->end($ORFstart-1); #modify start 
 				if ( ($strand == -1) or ($strand eq "-") ) {
 					$utr_feature->primary_tag('three_prime_utr');
 					create_or_replace_tag($utr_feature,'ID',$ID.'-utr3-'.$utr3_counter); #modify name
	 				push(@utr3_features, $utr_feature);#save that cds
	 				$utr3_counter++;
	 			}else{
	 				$utr_feature->primary_tag('five_prime_utr');
	 				create_or_replace_tag($utr_feature,'ID',$ID.'-utr5-'.$utr5_counter); #modify name
	 				push(@utr5_features, $utr_feature);#save that cds
	 				$utr5_counter++;
	 			}
 			}
 			if($exon_feature->end > $ORFend){
 				my $utr_feature=clone($exon_feature);#create a copy of the feature 
 				$utr_feature->start($ORFend+1); #modify start 
 				if ( ($strand == -1) or ($strand eq "-") ) {
 					$utr_feature->primary_tag('five_prime_utr');
	 				create_or_replace_tag($utr_feature,'ID',$ID.'-utr5-'.$utr5_counter); #modify name
	 				push(@utr5_features, $utr_feature);#save that cds
	 				$utr5_counter++;
	 			}else{
	 				$utr_feature->primary_tag('three_prime_utr');
 					create_or_replace_tag($utr_feature,'ID',$ID.'-utr3-'.$utr3_counter); #modify name
	 				push(@utr3_features, $utr_feature);#save that cds
	 				$utr3_counter++;
	 			}
 			}
		}
		# cds overlap fully an exon
		elsif( ($exon_feature->end <= $ORFend) and ($exon_feature->start >= $ORFstart) ){
 			my $cds_feature=clone($exon_feature);#create a copy of the feature 						exon    ========================
 			$cds_feature->primary_tag('cds');
 			#get old name 																			cds  ===============================
 			my @IDs = $cds_feature->get_tag_values('ID');
 			my $ID = $IDs[0];
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
 			$cds_feature->primary_tag('cds');
 				#get old name
 			my @IDs = $cds_feature->get_tag_values('ID');
 			my $ID = $IDs[0];
 			create_or_replace_tag($cds_feature,'ID',$ID.'-cds-'.$cds_counter); #modify name
 			push(@cds_features, $cds_feature);#save that cds
 			$cds_counter++;
 			#manage UTR
 			my $utr_feature=clone($exon_feature);#create a copy of the feature 
 			$utr_feature->start($ORFend+1); #modify end 
 			@IDs = $utr_feature->get_tag_values('ID');
	 		$ID = $IDs[0];
	 		if ( ($strand == -1) or ($strand eq "-") ) {
	 			$utr_feature->primary_tag('five_prime_utr');
	 			create_or_replace_tag($utr_feature,'ID',$ID.'-utr5-'.$utr5_counter); #modify name
	 			push(@utr5_features, $utr_feature);#save that cds
	 			$utr5_counter++;

	 		}else{
				$utr_feature->primary_tag('three_prime_utr');
	 			create_or_replace_tag($utr_feature,'ID',$ID.'-utr3-'.$utr3_counter); #modify name
	 			push(@utr3_features, $utr_feature);#save that cds
	 			$utr3_counter++;
	 		}

	      }
	      else{ #cds overlap start end exon
	       	#Manage CDS
	       	my $cds_feature=clone($exon_feature);#create a copy of the feature
 			$cds_feature->start($ORFstart); #modify start 										exon ===============================
 			$cds_feature->primary_tag('cds');	
 				#get old name 																					cds =====================================
 			my @IDs = $cds_feature->get_tag_values('ID');
 			my $ID = $IDs[0];
 			create_or_replace_tag($cds_feature,'ID',$ID.'-cds-'.$cds_counter); #modify name
 			push(@cds_features, $cds_feature);#save that cds
 			$cds_counter++;
	 		 #Manage UTR
 			my $utr_feature=clone($exon_feature);#create a copy of the feature 
 			$utr_feature->end($ORFstart-1); #modify start 
 			@IDs = $utr_feature->get_tag_values('ID');
	 		$ID = $IDs[0];
	 		if ( ($strand == -1) or ($strand eq "-") ) {
	 			$utr_feature->primary_tag('three_prime_utr');
	 			create_or_replace_tag($utr_feature,'ID',$ID.'-utr3-'.$utr3_counter); #modify name
	 			push(@utr3_features, $utr_feature);#save that cds
	 			$utr3_counter++;
	 		}else{
	 			$utr_feature->primary_tag('five_prime_utr');
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
	 			my @IDs = $utr_feature->get_tag_values('ID');
	 			my $ID = $IDs[0];
	 			if ( ($strand == -1) or ($strand eq "-") ) {
	 				$utr_feature->primary_tag('three_prime_utr');
	 				create_or_replace_tag($utr_feature,'ID',$ID.'-utr3-'.$utr3_counter); #modify name
	 				push(@utr3_features, $utr_feature);#save that cds
	 				$utr3_counter++;
	 			}else{
	 				$utr_feature->primary_tag('five_prime_utr');
	 				create_or_replace_tag($utr_feature,'ID',$ID.'-utr5-'.$utr5_counter); #modify name
	 				push(@utr5_features, $utr_feature);#save that cds
	 				$utr5_counter++;
	 			}
	 			
	 		}
	 		else{	#UTR3 in + strand
		    	my $utr_feature=clone($exon_feature);#create a copy of the feature 													exon ===============================
	 			#get old name
	 			my @IDs = $utr_feature->get_tag_values('ID'); 									#cds ===============================
	 			my $ID = $IDs[0];
	 			if ( ($strand == -1) or ($strand eq "-") ) {
	 				$utr_feature->primary_tag('five_prime_utr');
	 				create_or_replace_tag($utr_feature,'ID',$ID.'-utr5-'.$utr5_counter); #modify name
	 				push(@utr5_features, $utr_feature);#save that cds
	 				$utr5_counter++;
	 			}else{
	 				$utr_feature->primary_tag('three_prime_utr');
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
  if ($mRNA_feature->start != $exonStart){
    $mRNA_feature->start($exonStart);
  }
  elsif($mRNA_feature->end != $exonEnd){
    $mRNA_feature->end($exonEnd);
  }
}

# Check the start and end of gene feature based on its mRNAs;
sub check_gene_positions{

  my ($hash_omniscient, $gene_id)=@_;

  #####
  #Modify gene start-end (have to check size of each mRNA)
  my $geneExtremStart=1000000000000;
  my $geneExtremEnd=0;
  foreach my $primary_tag_key_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...

    foreach my $mrna_feature ( @{$hash_omniscient->{'level2'}{'mrna'}{lc($gene_id)}}) {
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
  my $gene_feature=$hash_omniscient->{'level1'}{'gene'}{lc($gene_id)};
  if ($gene_feature->start != $geneExtremStart){
      $gene_feature->start($geneExtremStart);
   }
  if($gene_feature->end != $geneExtremEnd){
    $gene_feature->end($geneExtremEnd);
  }
}

1;
