#!/usr/bin/perl -w

package BILS::GFF3::Statistics ;

use strict;
use warnings;
use Data::Dumper;
use Bio::SeqIO;
use Clone 'clone';
#use BILS::Handler::GXFhandler qw(:Ok);
use Exporter qw(import);

our $VERSION     = 1.00;
our @ISA         = qw(Exporter);
our @EXPORT_OK   = qw(gff3_statistics);
our %EXPORT_TAGS = ( DEFAULT => [qw()],
                 Ok    => [qw(gff3_statistics)]);
=head1 SYNOPSIS



=head1 DESCRIPTION

	A library to get statistics form gff3 file 
	
=cut	

# Calculate information necessary going through the omniscient only once
sub gff3_statistics {
	
	my ($hash_omniscient, $genome) = @_  ;
	
	my @result;

	#my $out = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 3);

	#check genome size

	my $genomeSize=undef;
	if($genome){
		if( $genome =~ /^[0-9]+$/){ #check if it's a number
			$genomeSize=$genome;
		}
		elsif($genome){
			my $seqio = Bio::SeqIO->new(-file => $genome, '-format' => 'Fasta');
			while(my $seq = $seqio->next_seq) {
	  			my $string = $seq->seq;
	  			$genomeSize += length($string);
	  		}
		}
	printf("%-45s%d%s", "Total sequence length", $genomeSize,"\n");
	}

	# get nb of each feature inomniscient;
	my %nb_feat;
	my %nb_spread_feat; #(utr and cds)
	my %size_feat;
	my %longest;
	my %shortest;
	my %tagList_relationship;
	my $utr_both_side=0;
	my %utr_at_least_one_side;

	foreach my $tag_l1 (keys %{$hash_omniscient->{'level1'}}){
		if (! $tagList_relationship{'level1'}{$tag_l1}){
			$tagList_relationship{'level1'}{$tag_l1}++;
		}

		foreach my $id (keys %{$hash_omniscient->{'level1'}{$tag_l1}}){
			my $feature_l1= $hash_omniscient->{'level1'}{$tag_l1}{$id};
			
			#count number of feature
			$nb_feat{$tag_l1}++;
			#compute feature size
			my $sizeFeature=($feature_l1->end-$feature_l1->start)+1;
	    	$size_feat{$tag_l1}+=$sizeFeature;
	    	# grab longest
	    	if ((! $longest{$tag_l1}) or ($longest{$tag_l1} < $sizeFeature)){
	    		$longest{$tag_l1}=$sizeFeature;
	    	}
	    	# grab shorter
	    	if ((! $shortest{$tag_l1}) or ($shortest{$tag_l1} > $sizeFeature)){
	    		$shortest{$tag_l1}=$sizeFeature;
	    	}

	    	foreach my $tag_l2 (keys %{$hash_omniscient->{'level2'}}){
	    		#append list of tag if necessary
	    		if (! exists ($tagList_relationship{'level2'}{$tag_l1}{$tag_l2})){
					$tagList_relationship{'level2'}{$tag_l1}{$tag_l2}++;
				}
	    		if(exists ($hash_omniscient->{'level2'}{$tag_l2}{$id})){
	    			foreach my $feature_l2 ( @{$hash_omniscient->{'level2'}{$tag_l2}{$id}} ){
	    				
	    				#count number of feature
	    				$nb_feat{$tag_l2}++;
						#compute feature size
	    				my $sizeFeature=($feature_l2->end-$feature_l2->start)+1;
	    	  			$size_feat{$tag_l2}+=$sizeFeature;
	    	  			# grab longest
	    	  			if ((! $longest{$tag_l2}) or ($longest{$tag_l2} < $sizeFeature)){
	    					$longest{$tag_l2}=$sizeFeature;
	    				}
	    				# grab shorter
	    				if ((! $shortest{$tag_l2}) or ($shortest{$tag_l2} > $sizeFeature)){
	    					$shortest{$tag_l2}=$sizeFeature;
	    				}

	    	  			my $id_l2=lc($feature_l2->_tag_value('ID'));
				    	foreach my $tag_l3 (keys %{$hash_omniscient->{'level3'}}){
				    		#append list of tag if necessary
				    		if (! exists( $tagList_relationship{'level3'}{$tag_l2}{$tag_l3})) {
								$tagList_relationship{'level3'}{$tag_l2}{$tag_l3}++;
							}

				    		if(exists ($hash_omniscient->{'level3'}{$tag_l3}{$id_l2})){
								my $sizeMultiFeat=0;
								my $feature;
				    			foreach my $feature_l3 ( @{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2}} ){
				    				$feature=$feature_l3;
				    				
				    				#compute feature size
				    				my $sizeFeature=($feature_l3->end-$feature_l3->start)+1;
				    	  			$size_feat{$tag_l3}+=$sizeFeature;
				    	  			
				    	  			if(($tag_l3 =~ /cds/) or ($tag_l3 =~ /utr/)){
				    	  				$sizeMultiFeat+=$sizeFeature;
				    	  				$nb_spread_feat{$tag_l3}++;
				    	  			}
				    	  			else{
				    	  				#count number of feature
				    					$nb_feat{$tag_l3}++;
					    	  			# grab longest
					    	  			if ((! $longest{$tag_l3}) or ($longest{$tag_l3} < $sizeFeature)){
		    								$longest{$tag_l3}=$sizeFeature;
		    							}
		    							# grab shorter
					    				if ((! $shortest{$tag_l3}) or ($shortest{$tag_l3} > $sizeFeature)){
					    					$shortest{$tag_l3}=$sizeFeature;
					    				}
					    			}

				    				#mange utr per mRNA
				    				if ($tag_l3 =~ /three_prime_utr/){
				    					if( exists ($hash_omniscient->{'level3'}{'five_prime_utr'}{$id_l2})){
				    						$utr_both_side++;
				    					}
				    				}
				    				if (($tag_l3 =~ /three_prime_utr/) or ($tag_l3 =~ /five_prime_utr/) ) {
				    					if (! exists ($utr_at_least_one_side{$id_l2}) ) {
				    						$utr_at_least_one_side{$id_l2}++;
				    					}
				    				}
				    			}

				    			#in that case the feature was split in several peaces that have been glue together
				    			if (($tag_l3 =~ /utr/) or ($tag_l3 =~ /cds/)){
				    				#count number of feature
				    				$nb_feat{$tag_l3}++;
				    				# grab longest
					    	  		if ((! $longest{$tag_l3}) or ($longest{$tag_l3} < $sizeMultiFeat)){
		    							$longest{$tag_l3}=$sizeMultiFeat;
		    						}
		    						# grab shorter
					    			if ((! $shortest{$tag_l3}) or ($shortest{$tag_l3} > $sizeMultiFeat)){
					    				$shortest{$tag_l3}=$sizeMultiFeat;
					    				#$out->write_feature($feature); 
					    			}
				    			}
				  			}
				  		}
				  	}
				}
			}
	    }
	}

	my $info_number = _info_number(\%tagList_relationship, \%nb_feat, $utr_both_side, \%utr_at_least_one_side);
	push @result, @$info_number;

	my $info_number_spread = info_number_spread(\%tagList_relationship, \%nb_spread_feat);
	push @result, @$info_number_spread;

	my $info_mean_per = _info_mean_per(\%tagList_relationship, \%nb_feat);
	push @result, @$info_mean_per;

	my $info_length = _info_length(\%tagList_relationship, \%size_feat);
	push @result, @$info_length;

	my $info_mean_length = _info_mean_length(\%tagList_relationship, \%nb_feat, \%size_feat);
	push @result, @$info_mean_length;

	if($genome){		
		my $info_coverage = _info_coverage(\%tagList_relationship, \%size_feat, $genomeSize);
		push @result, @$info_coverage;
	}
	
	my $info_longest = _info_longest(\%tagList_relationship, \%longest);
	push @result, @$info_longest;

	my $info_shortest = _info_shortest(\%tagList_relationship, \%shortest);
	push @result, @$info_shortest;

return \@result;
}

#####
# Give info about number of spread feature of each type
# spread feature can only be level3 (utr or cds)
sub info_number_spread{
	my ($tag_list, $nb_feat) = @_  ;
	
	my @resu;

	#print level3
	foreach my $tag_l2 (keys %{$tag_list->{'level3'}}){
		foreach my $tag_l3 (keys %{$tag_list->{'level3'}{$tag_l2}}){
			if($tag_l3 =~ /utr/ or $tag_l3 =~ /cds/){
		   		push @resu, sprintf("%-45s%d%s", "Number of exon of $tag_l3", $nb_feat->{$tag_l3},"\n");
		   	}
		}
	}
	return  \@resu;
}

#####
# Give info about number of feature of each type
sub _info_number {

	my ($tag_list, $nb_feat, $utr_both_side, $utr_at_least_one_side) = @_  ;
	my @resu;
	my $there_is_utr=undef;

	#print level1
	foreach my $tag_l1 (keys %{$tag_list->{'level1'}}){
		push @resu, sprintf("%-45s%d%s", "Number of $tag_l1"."s", $nb_feat->{$tag_l1},"\n");

		#print level2
	  	foreach my $tag_l2 (keys %{$tag_list->{'level2'}{$tag_l1}}){
		    push @resu, sprintf("%-45s%d%s", "Number of $tag_l2"."s", $nb_feat->{$tag_l2},"\n");

		    #print level3
		    foreach my $tag_l3 (keys %{$tag_list->{'level3'}{$tag_l2}}){
		    	if($tag_l3 =~ /utr/){$there_is_utr=1;}
		        push @resu, sprintf("%-45s%d%s", "Number of $tag_l3"."s", $nb_feat->{$tag_l3},"\n");
		    }
		}
	}
	#manage utr both side
	if($there_is_utr){
		push @resu, sprintf("%-45s%d%s", "Number of mrnas with utr both sides", $utr_both_side,"\n");
		my $nb_at_lest_one_side = keys %{$utr_at_least_one_side};
		push @resu, sprintf("%-45s%d%s", "Number of mrnas with at least one UTR", $nb_at_lest_one_side,"\n");
	}
	return \@resu;
}

#############
# Give info about shortest feature of each type
sub _info_shortest {

	my ($tag_list, $shortest) = @_  ;
	my @resu;

	#print level1
	foreach my $tag_l1 (keys %{$tag_list->{'level1'}}){
		push @resu, sprintf("%-45s%d%s", "Shortest $tag_l1"."s", $shortest->{$tag_l1},"\n");

		#print level2
	  	foreach my $tag_l2 (keys %{$tag_list->{'level2'}{$tag_l1}}){
		    push @resu, sprintf("%-45s%d%s", "Shortest $tag_l2"."s", $shortest->{$tag_l2},"\n");

		    #print level3
		    foreach my $tag_l3 (keys %{$tag_list->{'level3'}{$tag_l2}}){
		        push @resu, sprintf("%-45s%d%s", "Shortest $tag_l3"."s", $shortest->{$tag_l3},"\n");
		    }
		}
	}
	return \@resu;
}

#############
# Give info about longest feature of each type
sub _info_longest {

	my ($tag_list, $longest) = @_  ;
	my @resu;

	#print level1
	foreach my $tag_l1 (keys %{$tag_list->{'level1'}}){
		push @resu, sprintf("%-45s%d%s", "Longest $tag_l1"."s", $longest->{$tag_l1},"\n");

		#print level2
	  	foreach my $tag_l2 (keys %{$tag_list->{'level2'}{$tag_l1}}){
		    push @resu, sprintf("%-45s%d%s", "Longest $tag_l2"."s", $longest->{$tag_l2},"\n");

		    #print level3
		    foreach my $tag_l3 (keys %{$tag_list->{'level3'}{$tag_l2}}){
		       push @resu, sprintf("%-45s%d%s", "Longest $tag_l3"."s", $longest->{$tag_l3},"\n");
		    }
		}
	}
	return \@resu;
}

#############
# Give info about number mean of feature of each type per Parent type (mRNA per gene / cds per mRNA / etc)
sub _info_mean_per {
	my ($tag_list, $nb_feat) = @_  ;
	my @resu;

	#print level1
	foreach my $tag_l1 (keys %{$tag_list->{'level1'}}){

		#print level2
	  	foreach my $tag_l2 (keys %{$tag_list->{'level2'}{$tag_l1}}){
		    my $mean= $nb_feat->{$tag_l2}/$nb_feat->{$tag_l1};
	    	push @resu, sprintf("%-45s%.1f%s", "mean $tag_l2"."s per $tag_l1", $mean,"\n");

		    #print level3
		    foreach my $tag_l3 (keys %{$tag_list->{'level3'}{$tag_l2}}){
		       my $mean= $nb_feat->{$tag_l3}/$nb_feat->{$tag_l2};
	        	push @resu, sprintf("%-45s%.1f%s", "mean $tag_l3"."s per $tag_l2", $mean,"\n");
		    }
		}
	}
	return \@resu;
}
	
#############
# Give info about lenght of the total of features by type
sub _info_length {
	my ($tag_list, $size_feat) = @_  ;
	my @resu;

	#print level1
	foreach my $tag_l1 (keys %{$tag_list->{'level1'}}){
		push @resu, sprintf("%-45s%d%s", "Total $tag_l1 length", $size_feat->{$tag_l1},"\n");

		#print level2
	  	foreach my $tag_l2 (keys %{$tag_list->{'level2'}{$tag_l1}}){
		     push @resu, sprintf("%-45s%d%s", "Total $tag_l2 length", $size_feat->{$tag_l2},"\n");

		    #print level3
		    foreach my $tag_l3 (keys %{$tag_list->{'level3'}{$tag_l2}}){
		       push @resu, sprintf("%-45s%d%s", "Total $tag_l3 length", $size_feat->{$tag_l3},"\n");
		    }
		}
	}
	return \@resu;
}

#############
# Give info about mean lenght of features by type
sub _info_mean_length {
	my ($tag_list, $nb_feat, $size_feat) = @_  ;
	my @resu;

	#print level1
	foreach my $tag_l1 (keys %{$tag_list->{'level1'}}){
		my $meanl= $size_feat->{$tag_l1}/$nb_feat->{$tag_l1};
	    push @resu, sprintf("%-45s%d%s", "mean $tag_l1 length", $meanl,"\n");

		#print level2
	  	foreach my $tag_l2 (keys %{$tag_list->{'level2'}{$tag_l1}}){
		    my $meanl= $size_feat->{$tag_l2}/$nb_feat->{$tag_l2};
	    	push @resu, sprintf("%-45s%d%s", "mean $tag_l2 length", $meanl,"\n");

		    #print level3
		    foreach my $tag_l3 (keys %{$tag_list->{'level3'}{$tag_l2}}){
		    my $meanl= $size_feat->{$tag_l3}/$nb_feat->{$tag_l3};
	        push @resu, sprintf("%-45s%d%s", "mean $tag_l3 length", $meanl,"\n");
		    }
		}
	}
	return \@resu;
}

#############
# Give info about the features' coverage (by types) within/among the genome
sub _info_coverage {
	my ($tag_list,  $size_feat, $genomeSize) = @_  ;
	my @resu;

	#print level1
	foreach my $tag_l1 (keys %{$tag_list->{'level1'}}){
		my $perc= ($size_feat->{$tag_l1}*100)/$genomeSize;
	    push @resu, sprintf("%-45s%.1f%s", "% of genome covered by $tag_l1", $perc,"\n");

		#print level2
	  	foreach my $tag_l2 (keys %{$tag_list->{'level2'}{$tag_l1}}){
		    my $perc= ($size_feat->{$tag_l2}*100)/$genomeSize;
	    	push @resu, sprintf("%-45s%.1f%s", "% of genome covered by $tag_l2", $perc,"\n");

		    #print level3
		    foreach my $tag_l3 (keys %{$tag_list->{'level3'}{$tag_l2}}){
		    my $perc= ($size_feat->{$tag_l3}*100)/$genomeSize;
	        push @resu, sprintf("%-45s%.1f%s", "% of genome covered by $tag_l3", $perc,"\n");
		    }
		}
	}
	return \@resu;
}


1;
