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


	We create a complex hash of hash containing all information needeed.
	The data are scan from level 2 to level 1 and 3. We do that because different type of feature from level 2 can have same type of feature of level1. (e.g: mRNA => gene and tRNA => gene).
	So the structure of the hash created is the following: 
	{type_feature_level2}{'level'}{type_feature_level}{'flag'}='value';
	'level' can be level1, level2 or level accordingly, allow to go all over the data for printing by driving the data form level1 to level3.
	'flag' correspond to the type of information that has been saved in 'value'
	
=cut	

# Calculate information necessary going through the omniscient only once
# return a lisf of sub_list - Sub list contain all inforamtion level1,2,3 of all feature linked to a type of feature of level 2. 
# (eg: Gene(l1),mRNA(l2),cds(l3),exon(l3), where the type of level1 and level3 feature are only those linked to mRNA.)
sub gff3_statistics {
	
	my ($hash_omniscient, $genome) = @_  ;
	
	my @result_list;

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
	my %all_info;

	foreach my $tag_l2 (keys %{$hash_omniscient->{'level2'}}){
		foreach my $id_l1 (keys %{$hash_omniscient->{'level2'}{$tag_l2}}){
			my $one_f2 = $hash_omniscient->{'level2'}{$tag_l2}{$id_l1}[0];

			#######################
			#get feature1 and info
			my $feature_l1=undef;
			my $tag_l1;
			foreach my $tag_level1 (keys %{$hash_omniscient->{'level1'}}){
				if (exists ($hash_omniscient->{'level1'}{$tag_level1}{$id_l1})){
					$feature_l1=$hash_omniscient->{'level1'}{$tag_level1}{$id_l1};
					$tag_l1=$tag_level1;
					last;
				}
			}
			if(! $feature_l1){print "Problem ! We didnt retrieve the level1 feature with id $id_l1\n";exit;}
			#count number of feature
			$all_info{$tag_l2}{'level1'}{$tag_l1}{'nb_feat'}++;
			#compute feature size
			my $sizeFeature=($feature_l1->end-$feature_l1->start)+1;
			$all_info{$tag_l2}{'level1'}{$tag_l1}{'size_feat'}+=$sizeFeature;
	    	# grab longest
	    	if ((! $all_info{$tag_l2}{'level1'}{$tag_l1}{'longest'}) or ($all_info{$tag_l2}{'level1'}{$tag_l1}{'longest'} < $sizeFeature)){
	    		$all_info{$tag_l2}{'level1'}{$tag_l1}{'longest'}=$sizeFeature;
	    	}
	    	# grab shorter
	    	if ((! $all_info{$tag_l2}{'level1'}{$tag_l1}{'shortest'}) or ($all_info{$tag_l2}{'level1'}{$tag_l1}{'shortest'} > $sizeFeature)){
	    		$all_info{$tag_l2}{'level1'}{$tag_l1}{'shortest'}=$sizeFeature;
	    	}

	    	#####
	    	# get all level2
			foreach my $feature_l2 ( @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}} ){
				
				#count number of feature
				$all_info{$tag_l2}{'level2'}{$tag_l2}{'nb_feat'}++;
				#compute feature size
				my $sizeFeature=($feature_l2->end-$feature_l2->start)+1;
	  			$all_info{$tag_l2}{'level2'}{$tag_l2}{'size_feat'}+=$sizeFeature;
	  			# grab longest
	  			if ((! $all_info{$tag_l2}{'level2'}{$tag_l2}{'longest'}) or ($all_info{$tag_l2}{'level2'}{$tag_l2}{'longest'} < $sizeFeature)){
					$all_info{$tag_l2}{'level2'}{$tag_l2}{'longest'}=$sizeFeature;
				}
				# grab shorter
				if ((! $all_info{$tag_l2}{'level2'}{$tag_l2}{'shortest'}) or ($all_info{$tag_l2}{'level2'}{$tag_l2}{'shortest'} > $sizeFeature)){
					$all_info{$tag_l2}{'level2'}{$tag_l2}{'shortest'}=$sizeFeature;
				}

				######
				#get all level3
				my $utr3 = undef;
				my $utr5 = undef;
	  			my $id_l2=lc($feature_l2->_tag_value('ID'));
		    	foreach my $tag_l3 (keys %{$hash_omniscient->{'level3'}}){

		    		if(exists ($hash_omniscient->{'level3'}{$tag_l3}{$id_l2})){
						my $sizeMultiFeat=0;
						my $nb_multiFeature=0;
						#Initialize intron to 0 to avoid error during printing results
						if(! exists ($all_info{$tag_l2}{'level3'}{$tag_l3}{'intron'})){$all_info{$tag_l2}{'level3'}{$tag_l3}{'intron'}=0;}

		    			foreach my $feature_l3 ( @{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2}} ){
		    				
		    				#count number feature of tag_l3 type
		    				$nb_multiFeature++;

		    				#compute feature size
		    				my $sizeFeature=($feature_l3->end-$feature_l3->start)+1;
		    				$all_info{$tag_l2}{'level3'}{$tag_l3}{'size_feat'}+=$sizeFeature;
		    	  			
		    	  			if(($tag_l3 =~ /cds/) or ($tag_l3 =~ /utr/)){
		    	  				$sizeMultiFeat+=$sizeFeature;
		    	  				$all_info{$tag_l2}{'level3'}{$tag_l3}{'nb_spread_feat'}++;
		    	  			}
		    	  			else{
		    	  				#count number of feature
		    					$all_info{$tag_l2}{'level3'}{$tag_l3}{'nb_feat'}++;
			    	  			# grab longest
			    	  			if ((! $all_info{$tag_l2}{'level3'}{$tag_l3}{'longest'}) or ($all_info{$tag_l2}{'level3'}{$tag_l3}{'longest'} < $sizeFeature)){
    								$all_info{$tag_l2}{'level3'}{$tag_l3}{'longest'}=$sizeFeature;
    							}
    							# grab shorter
			    				if ((! $all_info{$tag_l2}{'level3'}{$tag_l3}{'shortest'}) or ($all_info{$tag_l2}{'level3'}{$tag_l3}{'shortest'} > $sizeFeature)){
			    					$all_info{$tag_l2}{'level3'}{$tag_l3}{'shortest'}=$sizeFeature;
			    				}
			    			}
			    			####################
		    				#mange utr per mRNA
		    				if ($tag_l3 =~ /three_prime_utr/){
								$utr3=1;
		    				}
		    				if ($tag_l3 =~ /five_prime_utr/){
		    					$utr5=1;
		    				}
		    			}

		    			#in that case the feature was split in several peaces that have been glue together
		    			if (($tag_l3 =~ /utr/) or ($tag_l3 =~ /cds/)){
		    				#count number of feature
		    				$all_info{$tag_l2}{'level3'}{$tag_l3}{'nb_feat'}++;
		    				# grab longest
			    	  		if ((! $all_info{$tag_l2}{'level3'}{$tag_l3}{'longest'}) or ($all_info{$tag_l2}{'level3'}{$tag_l3}{'longest'} < $sizeMultiFeat)){
    							$all_info{$tag_l2}{'level3'}{$tag_l3}{'longest'}=$sizeMultiFeat;
    						}
    						# grab shorter
			    			if ((! $all_info{$tag_l2}{'level3'}{$tag_l3}{'shortest'}) or ($all_info{$tag_l2}{'level3'}{$tag_l3}{'shortest'} > $sizeMultiFeat)){
			    				$all_info{$tag_l2}{'level3'}{$tag_l3}{'shortest'}=$sizeMultiFeat;
			    			}
		    			}

		    			#Manage number intron per type of tag_l3
		    			if($nb_multiFeature > 1){
		    				my $nbIntron = ($nb_multiFeature-1);
		    				$all_info{$tag_l2}{'level3'}{$tag_l3}{'intron'}+=$nbIntron;
		    			}
		  			}
		  		}# END all feature level 3

		    	# 1) Manage UTR both side 
		    	if ($utr3  and $utr5){
			    		$all_info{$tag_l2}{'level2'}{$tag_l2}{'utr_both_side'}++;
		    	} # 2) Manage UTR at least one side 
		    	elsif ($utr3  or $utr5){ 
		   				$all_info{$tag_l2}{'level2'}{$tag_l2}{'utr_at_least_one_side'}++;
 				}
		  	}
		}
	}

	foreach my $type (keys %all_info){
		my $hashType = $all_info{$type};
		my @result;

		my $info_number = _info_number($hashType);
		push @result, @$info_number;

		my $info_number_spread = info_number_spread($hashType);
		push @result, @$info_number_spread;

		my $info_mean_per = _info_mean_per($hashType);
		push @result, @$info_mean_per;

		my $info_length = _info_length($hashType);
		push @result, @$info_length;

		my $info_mean_length = _info_mean_length($hashType);
		push @result, @$info_mean_length;

		if($genome){		
			my $info_coverage = _info_coverage($hashType, $genomeSize);
	 		push @result, @$info_coverage;
		}
	
		my $info_longest = _info_longest($hashType);
		push @result, @$info_longest;

		my $info_shortest = _info_shortest($hashType);
		push @result, @$info_shortest;

		push @result_list, \@result;
	}

return \@result_list;
}

#####
# Give info about number of spread feature of each type
# spread feature can only be level3 (utr or cds)
sub info_number_spread{
	my ($all_info) = @_  ;
	
	my @resu;

	#print level3
	foreach my $tag_l3 (keys %{$all_info->{'level3'}}){
	    #manage nb_spread_feat
	    if(exists ($all_info->{'level3'}{$tag_l3}{'nb_spread_feat'})){
	    	push @resu, sprintf("%-45s%d%s", "Number of exon of $tag_l3", $all_info->{'level3'}{$tag_l3}{'nb_spread_feat'},"\n");
	    	push @resu, sprintf("%-45s%d%s", "Number of intron of $tag_l3", $all_info->{'level3'}{$tag_l3}{'intron'},"\n");
	    }
	}
	return  \@resu;
}

#####
# Give info about number of feature of each type
sub _info_number {

	my ($all_info) = @_  ;
	my @resu;
	my $there_is_utr=undef;

	#print level1
	foreach my $tag_l1 (keys %{$all_info->{'level1'}}){
		push @resu, sprintf("%-45s%d%s", "Number of $tag_l1"."s", $all_info->{'level1'}{$tag_l1}{'nb_feat'},"\n");
	}

	#print level2
	foreach my $tag_l2 (keys %{$all_info->{'level2'}}){
	    push @resu, sprintf("%-45s%d%s", "Number of $tag_l2"."s", $all_info->{'level2'}{$tag_l2}{'nb_feat'},"\n");
	    #manage utr both side
	    if(exists ($all_info->{'level2'}{$tag_l2}{'utr_both_side'})){
			push @resu, sprintf("%-45s%d%s", "Number of mrnas with utr both sides", $all_info->{'level2'}{$tag_l2}{'utr_both_side'},"\n");
		}
		#manage utr both side
		if(exists ($all_info->{'level2'}{$tag_l2}{'utr_at_least_one_side'})){
			push @resu, sprintf("%-45s%d%s", "Number of mrnas with at least one utr", $all_info->{'level2'}{$tag_l2}{'utr_at_least_one_side'},"\n");
		}
	 }

	#print level3
	foreach my $tag_l3 (keys %{$all_info->{'level3'}}){
	    push @resu, sprintf("%-45s%d%s", "Number of $tag_l3"."s", $all_info->{'level3'}{$tag_l3}{'nb_feat'},"\n");
	    #intron case 
	    if($tag_l3 eq "exon"){
	    	push @resu, sprintf("%-45s%d%s", "Number of intron"."s", $all_info->{'level3'}{$tag_l3}{'intron'},"\n");
	    }	    
	}
	

	return \@resu;
}

#############
# Give info about shortest feature of each type
sub _info_shortest {

	my ($all_info) = @_  ;
	my @resu;

	#print level1
	foreach my $tag_l1 (keys %{$all_info->{'level1'}}){
		push @resu, sprintf("%-45s%d%s", "Shortest $tag_l1"."s", $all_info->{'level1'}{$tag_l1}{'shortest'},"\n");
	}

	#print level2
	foreach my $tag_l2 (keys %{$all_info->{'level2'}}){
	    push @resu, sprintf("%-45s%d%s", "Shortest $tag_l2"."s", $all_info->{'level2'}{$tag_l2}{'shortest'},"\n");
	 }

	#print level3
	foreach my $tag_l3 (keys %{$all_info->{'level3'}}){
	    push @resu, sprintf("%-45s%d%s", "Shortest $tag_l3"."s", $all_info->{'level3'}{$tag_l3}{'shortest'},"\n");
	}

	return \@resu;
}

#############
# Give info about longest feature of each type
sub _info_longest {

	my ($all_info) = @_  ;
	my @resu;

	#print level1
	foreach my $tag_l1 (keys %{$all_info->{'level1'}}){
		push @resu, sprintf("%-45s%d%s", "Longest $tag_l1"."s", $all_info->{'level1'}{$tag_l1}{'longest'},"\n");
	}

	#print level2
	foreach my $tag_l2 (keys %{$all_info->{'level2'}}){
	    push @resu, sprintf("%-45s%d%s", "Longest $tag_l2"."s", $all_info->{'level2'}{$tag_l2}{'longest'},"\n");
	 }

	#print level3
	foreach my $tag_l3 (keys %{$all_info->{'level3'}}){
	    push @resu, sprintf("%-45s%d%s", "Longest $tag_l3"."s", $all_info->{'level3'}{$tag_l3}{'longest'},"\n");
	}

	return \@resu;
}

#############
# Give info about number mean of feature of each type per Parent type (mRNA per gene / cds per mRNA / etc)
sub _info_mean_per {
	my ($all_info) = @_  ;
	my @resu;

	#print level2
	foreach my $tag_l2 (keys %{$all_info->{'level2'}}){
		foreach my $tag_l1 (keys %{$all_info->{'level1'}}){
			my $mean=  $all_info->{'level2'}{$tag_l2}{'nb_feat'}/$all_info->{'level1'}{$tag_l1}{'nb_feat'};
		    push @resu, sprintf("%-45s%.1f%s", "mean $tag_l2"."s per $tag_l1", $mean,"\n");
		}
	}
	#print level3
	foreach my $tag_l3 (keys %{$all_info->{'level3'}}){
	    foreach my $tag_l2 (keys %{$all_info->{'level2'}}){
			my $mean=  $all_info->{'level3'}{$tag_l3}{'nb_feat'}/$all_info->{'level2'}{$tag_l2}{'nb_feat'};
		    push @resu, sprintf("%-45s%.1f%s", "mean $tag_l3"."s per $tag_l2", $mean,"\n");
		}
	}

	return \@resu;
}
	
#############
# Give info about lenght of the total of features by type
sub _info_length {
	my ($all_info) = @_  ;
	my @resu;

	#print level1
	foreach my $tag_l1 (keys %{$all_info->{'level1'}}){
		push @resu, sprintf("%-45s%d%s", "Total $tag_l1 length", $all_info->{'level1'}{$tag_l1}{'size_feat'},"\n");
	}

	#print level2
	foreach my $tag_l2 (keys %{$all_info->{'level2'}}){
	    push @resu, sprintf("%-45s%d%s", "Total $tag_l2 length", $all_info->{'level2'}{$tag_l2}{'size_feat'},"\n");
	 }

	#print level3
	foreach my $tag_l3 (keys %{$all_info->{'level3'}}){
	    push @resu, sprintf("%-45s%d%s", "Total $tag_l3 length", $all_info->{'level3'}{$tag_l3}{'size_feat'},"\n");
	}

	return \@resu;
}

#############
# Give info about mean lenght of features by type
sub _info_mean_length {
	my ($all_info) = @_  ;
	my @resu;


	#print level1
	foreach my $tag_l1 (keys %{$all_info->{'level1'}}){
		my $meanl= $all_info->{'level1'}{$tag_l1}{'size_feat'}/$all_info->{'level1'}{$tag_l1}{'nb_feat'};
		push @resu, sprintf("%-45s%d%s", "mean $tag_l1 length", $meanl,"\n");
	}

	#print level2
	foreach my $tag_l2 (keys %{$all_info->{'level2'}}){
		my $meanl= $all_info->{'level2'}{$tag_l2}{'size_feat'}/$all_info->{'level2'}{$tag_l2}{'nb_feat'};
	    push @resu, sprintf("%-45s%d%s", "mean $tag_l2 length", $meanl,"\n");
	 }

	#print level3
	foreach my $tag_l3 (keys %{$all_info->{'level3'}}){
		my $meanl= $all_info->{'level3'}{$tag_l3}{'size_feat'}/$all_info->{'level3'}{$tag_l3}{'nb_feat'};
	    push @resu, sprintf("%-45s%d%s", "mean $tag_l3 length", $meanl,"\n");
	}

	return \@resu;
}

#############
# Give info about the features' coverage (by types) within/among the genome
sub _info_coverage {
	my ($all_info, $genomeSize) = @_  ;
	my @resu;

	#print level1
	foreach my $tag_l1 (keys %{$all_info->{'level1'}}){
		my $perc= ($all_info->{'level1'}{$tag_l1}{'size_feat'}*100)/$genomeSize;
		push @resu, sprintf("%-45s%.1f%s", "% of genome covered by $tag_l1", $perc,"\n");
	}

	#print level2
	foreach my $tag_l2 (keys %{$all_info->{'level2'}}){
		my $perc= ($all_info->{'level2'}{$tag_l2}{'size_feat'}*100)/$genomeSize;
	    push @resu, sprintf("%-45s%.1f%s", "% of genome covered by $tag_l2", $perc,"\n");
	 }

	#print level3
	foreach my $tag_l3 (keys %{$all_info->{'level3'}}){
	    my $perc= ($all_info->{'level3'}{$tag_l3}{'size_feat'}*100)/$genomeSize;
	    push @resu, sprintf("%-45s%.1f%s", "% of genome covered by $tag_l3", $perc,"\n");
	}

	return \@resu;
}


1;
