#!/usr/bin/perl -w

package BILS::Handler::GXFhandler ;

use strict;
use warnings;
use Data::Dumper;
use File::Basename;
use Try::Tiny;
use LWP::UserAgent;
use Bio::OntologyIO::obo;
use Clone 'clone';
use BILS::Handler::GFF3handler qw(:Ok);
use Exporter qw(import);
use Bio::Tools::GFF;

our $VERSION     = 1.00;
our @ISA         = qw(Exporter);
our @EXPORT_OK   = qw(get_level select_gff_format check_mrna_positions modelate_utr_and_cds_features_from_exon_features_and_cds_start_stop slurp_gff3_file_JD _check_all_level1_positions _check_all_level2_positions);
our %EXPORT_TAGS = ( DEFAULT => [qw()],
                 	 Ok    => [qw(get_level select_gff_format check_mrna_positions modelate_utr_and_cds_features_from_exon_features_and_cds_start_stop slurp_gff3_file_JD _check_all_level1_positions _check_all_level2_positions)]);
=head1 SYNOPSIS


=head1 DESCRIPTION

	A library to convert handle any kind of gff file and save it in memory as GFF3 "full" format.
	Full format means, we expand exon with several parent, we add ID everywhere (even if level3 ID is not mandatory), and Parent everywhere.

=head1 VERSION
  
    Perl librairy last edited 2-Nov-2016.

=head1 CONTACT
    jacques.dainat@bils.se (Jacques Dainat)

=cut	

#===== TO  DO ===== 
# When creating a parent check its type from the value in the constant hash

#===== INFO ===== 	
# _manage_location and _check_overlap_name_diff are methods which compare list of list to look at overlap* (next ot each other is considered as overlap for the first method - mandatory). The first method
# rechek all the element if one modification is done, the second one go through element from left to right.  


##########################
#     DEFINE CONSTANT    # 
#/!\ Has to be defined in lowercase !
#level1 are features without parent
use constant LEVEL1 => { gene => 1, ncrna_gene => 2, mirna_gene => 3, lincrna_gene => 4, rrna_gene => 5, snrna_gene => 6, snorna_gene => 7, pseudogene => 8, sirna_gene => 10, pirna_gene => 11, sts => 12,
						 match => 9, protein_match => 13,  cDNA_match => 14, EST_match => 15, translated_nucleotide_match => 16, nucleotide_to_protein_match => 17, nucleotide_motif => 18, polypeptide => 19};
#level2 are features with parent and potentially a child (no child as example for match)
use constant LEVEL2 => { mrna => "gene", transcript => "gene", processed_transcript => "gene", nmd_transcript_variant => "gene", aberrant_processed_transcript => "gene", rna => "gene", trna => "gene",
						ncrna => "ncrna_gene", mirna => "mirna_gene", lcrna => "lincrna_gene", lincrna => "lincrna_gene", rrna => "rrna_gene", snrna => "snrna_gene", snorna => "snorna_gene",  
						trna_pseudogene => "pseudogene", pseudogenic_transcript => "pseudogene", processed_pseudogene => "pseudogene", unitary_pseudogene => "pseudogene",
						match_part => "match", similarity => "match", 
						nc_primary_transcript => "gene",
						misc_rna => "gene", rnase_p_rna => "gene", tmrna => "gene", srp_rna => "gene", vaultrna => "gene", sirna => "sirna_gene", piRNA => "pirna_gene"};
#level3 features have parents, but no child
# cds => "exon" mean that cds is included into an exon				
use constant LEVEL3 => { cds => "exon", exon => 2, stop_codon => "exon", start_codon => "exon", three_prime_utr => "exon", five_prime_utr => "exon", utr => "exon", selenocysteine => 8, non_canonical_three_prime_splice_site => 8, non_canonical_five_prime_splice_site => 10,
						stop_codon_read_through => "exon", sig_peptide => "exon", tss => "exon", tts => "exon", intron => 15 };
#feature that can be split over different locations
use constant SPREADFEATURE => {cds => 1, three_prime_utr => 2, five_prime_utr => 3, utr => 4};
use constant PREFIXL2 => "nbis_noL2id";

# ====== PURPOSE =======:
# Save in omniscient hash (sorted in a specific way (3 levels)) a whole gff3 file
# Parser phylosophy: Parse by Parent/child ELSE 
#						Parse by comon_tag  ELSE
#							Parse by sequential (mean group feattures in a bucket, and the bucket change at each level2 feature, and bucket are join in a comon tag at each new L1 feature)
#Priority Parent > locus_tag > sequential
# ====== INPUT =======:
# $file => string (file)
# $comonTagAttribute => list of tags to consider for gathering features
# $gffVersion => Int (if is used, force the parser to use this gff parser instead of guessing)
# $verbose =>define the deepth of verbosity
# $nocheck is to deactivate sanity check. It's in devellopement. We should be able to tune the deactivation of selected checks.
my $createL3forL2orphan = 1;
my $fh_error = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => 3);
sub slurp_gff3_file_JD {
	
	my ($self, $file, $comonTagAttribute, $gffVersion, $verbose, $nocheck, $quiet) = @_  ;

	## TO DO 
	# Handle the parameter with a hash
	if(! $verbose){$verbose=0;}
	my $overlapCheck = 0; # option to activate the check for overlapping locus

	my $start_run = time();
	my $previous_time = undef;	

	#get header information 
	my $gff3headerInfo = _check_header($file);

#	+-----------------------------------------+
#	|	HANDLE SOFA (feature-ontology)		  |
#	+-----------------------------------------+

	my $ontology = {};
	my $ontology_obj = _handle_ontology($gff3headerInfo, $verbose, $quiet);
	if($ontology_obj){
		$ontology = create_term_and_id_hash($ontology_obj);
	}

	if(! keys %{$ontology} ){ #hash is empty
		print "No data retrieved among the feature-ontology.\n" unless $quiet;
	}

#	+-----------------------------------------+
#	|			HANDLE WARNING 				  |
#	+-----------------------------------------+
	# Handle to not print to much warning
	my %WARNS;
	my %globalWARNS;
	my $nbWarnLimit=100;
  	local $SIG{__WARN__} = sub {
    my $message = shift;
    my @thematic=split /@/,$message ;
    
    if($thematic[0] eq "GLOBAL"){ #extract global warning
    	push @{$globalWARNS{$thematic[1]}}, $thematic[2];
    }
    else{
	    $WARNS{$thematic[0]}++;
		    if($verbose >= 1){
		    	if ($WARNS{$thematic[0]} <= $nbWarnLimit){
					print $message;
				}
				if($WARNS{$thematic[0]} == $nbWarnLimit){
					print "$thematic[0] ************** Too much WARNING message we skip the next **************\n";
				}
		  	}
	  	}
	 };

#	+-----------------------------------------+
#	|			HANDLE FEATUTRES PARSING ACCORDING TO TYPE 			  |
#	+-----------------------------------------+

	my %mRNAGeneLink; #Hast that keep track about link between l2 and l1
	my %omniscient; #Hast where all the feature will be saved
	my %duplicate;# Hash to store duplicated feature info
	my %miscCount;# Hash to store any counter. Will be use to create a new uniq ID
	my %uniqID;# Hash to follow up with an uniq identifier every feature
	my %uniqIDtoType; # Hash to follow up with an uniq identifier every feature type
	my %locusTAG;
	my %infoSequential;# Hash to store
	my $locusTAGvalue=undef;
	my $last_l1_f=undef;
	my $last_l2_f=undef;
	my $last_l3_f=undef;
	my $last_f=undef;# last feature handled
	my $lastL1_new =undef; # Bolean to check if last l1 feature is a newly created one. Important to deal with strict sequential

	# ARRAY CASE <============================
	if(ref($file) eq 'ARRAY'){
		 foreach my $feature (@{$file}) {
			($locusTAGvalue, $last_l1_f, $last_l2_f, $last_l3_f, $last_f, $lastL1_new) = manage_one_feature($ontology, $feature, \%omniscient, \%mRNAGeneLink, \%duplicate, \%miscCount, \%uniqID, \%uniqIDtoType, \%locusTAG, \%infoSequential, $comonTagAttribute, $locusTAGvalue, $last_l1_f, $last_l2_f, $last_l3_f, $last_f, $lastL1_new, $verbose);
   		}
	}
	# HASH CASE <============================
	elsif(ref($file) eq 'HASH'){

		foreach my $level (keys %{$file}){
			if ( ref($file->{$level}) eq 'HASH'){ #level1,level2,#level3
				foreach my $tag (keys %{$file->{$level}}){ 
					foreach my $id (keys %{$file->{$level}{$tag}}){ 
						if ( ref($file->{$level}{$tag}{$id}) eq 'ARRAY'){ #level2,#level3
							foreach my $feature ( @{$file->{$level}{$tag}{$id} }){ 
								($locusTAGvalue, $last_l1_f, $last_l2_f, $last_l3_f, $last_f, $lastL1_new) = manage_one_feature($ontology, $feature, \%omniscient, \%mRNAGeneLink, \%duplicate, \%miscCount, \%uniqID, \%uniqIDtoType, \%locusTAG, \%infoSequential, $comonTagAttribute, $locusTAGvalue, $last_l1_f, $last_l2_f, $last_l3_f, $last_f, $lastL1_new, $verbose);
							}
						}
						else{ #level1
							my $feature = $file->{$level}{$tag}{$id};
							($locusTAGvalue, $last_l1_f, $last_l2_f, $last_l3_f, $last_f, $lastL1_new) = manage_one_feature($ontology, $feature, \%omniscient, \%mRNAGeneLink, \%duplicate, \%miscCount, \%uniqID, \%uniqIDtoType, \%locusTAG, \%infoSequential, $comonTagAttribute, $locusTAGvalue, $last_l1_f, $last_l2_f, $last_l3_f, $last_f, $lastL1_new, $verbose);
						}
					}
				}
			}
			else{ #extra list of feature
				foreach my $feature ( @{$file->{$level}} ) {
					($locusTAGvalue, $last_l1_f, $last_l2_f, $last_l3_f, $last_f, $lastL1_new) = manage_one_feature($ontology, $feature, \%omniscient, \%mRNAGeneLink, \%duplicate, \%miscCount, \%uniqID, \%uniqIDtoType, \%locusTAG, \%infoSequential, $comonTagAttribute, $locusTAGvalue, $last_l1_f, $last_l2_f, $last_l3_f, $last_f, $lastL1_new, $verbose);
   				}
			}
		}
	}
	# FILE CASE <============================
	else{

		#GFF format used for parser
		my $format;
		if($gffVersion){$format = $gffVersion;}
		else{ $format = select_gff_format($file);}

		print "=>GFF version parser used: $format\n";
		my $gffio = Bio::Tools::GFF->new(-file => $file, -gff_version => $format);
		
		#read every lines
		while( my $feature = $gffio->next_feature()) {
			if($format eq "1"){_gff1_corrector($feature);} # case where gff1 has been used to parse.... we have to do some attribute manipulations
			($locusTAGvalue, $last_l1_f, $last_l2_f, $last_l3_f, $last_f, $lastL1_new) = manage_one_feature($ontology, $feature, \%omniscient, \%mRNAGeneLink, \%duplicate, \%miscCount, \%uniqID, \%uniqIDtoType, \%locusTAG, \%infoSequential, $comonTagAttribute, $locusTAGvalue, $last_l1_f, $last_l2_f, $last_l3_f, $last_f, $lastL1_new, $verbose);
	    }

	    #close the file
	    $gffio->close();
	}

    if($verbose  >= 1) {_printSurrounded("parsing done in ".(time() - $start_run)." seconds",30,".","\n"); $previous_time = time();}

#	+-----------------------------------------+
#	|			CHECK OMNISCIENT			  |
#	+-----------------------------------------+

    _printSurrounded("Start extra check",50,"#","\n") if $verbose;
    _printSurrounded("Check0: _check_duplicates",30,"*","\n") if ($verbose >= 1) ;

    #Check if duplicate detected:
    _check_duplicates(\%duplicate, \%omniscient, $verbose);
	_printSurrounded("Check1: _check_sequential",30,"*") if ($verbose >= 1) ;

print Dumper(\%omniscient);

	if(! $nocheck){
	    #Check sequential if we can fix cases. Hash to be done first, else is risky that we remove orphan L1 feature ... that are not yet linked to a sequential bucket
		if( keys %infoSequential ){ #hash is not empty
	    	_check_sequential(\%infoSequential, \%omniscient, \%miscCount, \%uniqID, \%uniqIDtoType, \%locusTAG, \%mRNAGeneLink, $verbose);
	    	undef %infoSequential;
	    }
	
	    else{ print "Nothing to check as sequential !\n" if($verbose >= 1) }
		if($verbose >= 1) {print "      done in ",time() - $previous_time," seconds\n\n\n"; $previous_time = time();}
		_printSurrounded("Check2: _remove_orphan_l1",30,"*") if($verbose >= 1) ;
	}

	    #check level1 has subfeature else we remove it
	  	_remove_orphan_l1(\%omniscient, \%miscCount, \%uniqID, \%uniqIDtoType, \%mRNAGeneLink, $verbose); #or fix if level2 is missing (refseq case)
		if($verbose >= 1) {print "      done in ",time() - $previous_time," seconds\n\n\n" ; $previous_time = time();}
	  	_printSurrounded("Check3: _check_l3_link_to_l2",30,"*") if ($verbose >= 1) ;

	if(! $nocheck){
	    #Check relationship between l3 and l2
	    _check_l3_link_to_l2(\%omniscient, \%mRNAGeneLink, \%miscCount, \%uniqID, \%uniqIDtoType, $verbose); # When creating L2 missing we create as well L1 if missing too
		if($verbose >= 1) {print "      done in ",time() - $previous_time," seconds\n\n\n" ; $previous_time = time();}
		_printSurrounded("Check4: _check_exons",30,"*") if ($verbose >= 1) ;

	    #Check relationship L3 feature, exons have to be defined... / mRNA position are checked!
	    _check_exons(\%omniscient, \%mRNAGeneLink, \%miscCount, \%uniqID,  \%uniqIDtoType, $verbose);
		if($verbose >= 1) {print "      done in ",time() - $previous_time," seconds\n\n\n"; $previous_time = time();}
		_printSurrounded("Check5: _check_utrs",30,"*") if ($verbose >= 1) ;

		#Check relationship L3 feature, exons have to be defined... / mRNA position are checked!
	    _check_utrs(\%omniscient, \%mRNAGeneLink, \%miscCount, \%uniqID,  \%uniqIDtoType, $verbose);
		if($verbose >= 1) {print "      done in ",time() - $previous_time," seconds\n\n\n"; $previous_time = time();}
		_printSurrounded("Check6: _check_gene_link_to_mrna",30,"*") if ($verbose >= 1) ;

	    #Check relationship between mRNA and gene.  / gene position are checked! If No Level1 we create it !
	    _check_gene_link_to_mrna(\%omniscient, $verbose);
		if($verbose >= 1) {print "      done in ",time() - $previous_time," seconds\n\n\n" ; $previous_time = time();}
		_printSurrounded("Check7: _check_all_level2_positions",30,"*") if ($verbose >= 1) ;

		# Check rna positions compared to its l2 features
		_check_all_level2_positions(\%omniscient, $verbose);
		if($verbose >= 1) {print "      done in ",time() - $previous_time," seconds\n\n\n" ; $previous_time = time();}
		_printSurrounded("Check8: _check_all_level1_positions",30,"*") if ($verbose >= 1) ;

		# Check gene positions compared to its l2 features
		_check_all_level1_positions(\%omniscient, $verbose);
		if($verbose >= 1) {print "      done in ",time() - $previous_time," seconds\n\n\n" ; $previous_time = time();}
		

		#check loci names (when overlap should be the same if type is the same)
		if ($overlapCheck){
			_printSurrounded("Check9: _check_overlap_name_diff",30,"*") if ($verbose >= 1) ;
			_check_overlap_name_diff(\%omniscient, \%mRNAGeneLink, $verbose);
		    if($verbose >= 1)  {print "      done in ",time() - $previous_time," seconds\n\n\n" ; $previous_time = time();}
		}

		#check identical isoforms
		_printSurrounded("Check10: _check_identical_isoforms",30,"*") if ($verbose >= 1) ;
		_check_identical_isoforms(\%omniscient, \%mRNAGeneLink, $verbose);
		 if($verbose >= 1)  {print "      done in ",time() - $previous_time," seconds\n\n\n" ; $previous_time = time();}
	}

    # To keep track of How many Warnings we got....
    foreach my $thematic (keys %WARNS){
  		my $nbW = $WARNS{$thematic};
  		if($nbW > $nbWarnLimit){
  			print "$nbW warning messages: $thematic\n";
  		}	
  	}

	print "Parsing and check done in ", time() - $start_run," seconds\n\n\n" if ($verbose >= 1);

#	+-----------------------------------------+
#	|			HANDLE GLOBAL WARNING 		  |
#	+-----------------------------------------+
	if( keys %globalWARNS ){
		if(exists($globalWARNS{"parser1"}) ) {
			my %hash   = map { $_, 1 } @{$globalWARNS{parser1}};
			my @unique = keys %hash;
			my $string = "Primary tag values (3rd column) not expected => @unique\n".
			"Those primary tag are not yet taken into account !\n".
			"Please modify the code to define one of the three possible levels of this feature or provide this information on the fly through the corresponding parameter (level1, level2 or level3).\n".
			"To resume:\n".
			"- it must be a level1 feature if it has no parent.\n".
			"- it must be a level2 feature if it has a parent and this parent is from level1.\n".
			"- it must be a level3 feature if it has a parent and this parent has also a parent.\n\n".
			"Currently the tool just ignore them, So if they where Level1,level2, a gene or RNA feature will be created accordingly.";
			_printSurrounded($string,150,"*") ;
		}
		if(exists($globalWARNS{"ontology1"}) ) {
			if( keys %{$ontology} ){
				my $string = "Primary tag values (3rd column) not expected => @{$globalWARNS{ontology1}}\n".
				"Those values are not compatible with gff3 format and the tool cannot guess to which one they correspond to.\n".
				"If you want to follow rigourously the gff3 format, please visit this website:\n".
				"https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md\n".
				"They provide tools to check the gff3 format.\n".
				"Even if you have this warning, you should be able to use the gff3 output in most of gff3 tools.";
				_printSurrounded($string,150,"*") ;
			}
			else{
				my $string = "No feature-ontology was available, we haven't checked that the feature types (3rd column) correspond to the gff3 specifactions.";
				_printSurrounded($string,150,"*") ;
			}
		}
	}

    #return
	return \%omniscient, \%mRNAGeneLink	;
}

# ====== PURPOSE =======:
# The method read a gff3 feature, Check for the sanity according to what will has been read before, and the whole features that will be read.
# Designed to be used within a loop going through a big amout of feature
# ====== INPUT =======:
# $feature => gff feature object
# 		example: scaffold625	maker	exon	341518	341628	.	+	.	ID=CLUHART00000008717:exon:1406;Parent=CLUHART00000008717
# $omniscient => hash to store all the gff feature in 3 levels structures
# 		example at level1: $omniscient->{"level1"}{$primary_tag}{$id}=$feature;
# 		example at other levels: $omniscient->{"levelX"}{$primary_tag}{$parent}= [$feature];
# $mRNAGeneLink => hash to keep track link from L2 to l1 (avoid to go through the whole omniscient to retrieve this information)
# 		example: $mRNAGeneLink->{lc($id)}=$parent;
# $duplicate => hash to keep track of duplicates found
# 		example: duplicate->{$level}{$primary_tag}{$id} = [$feature];
# $miscCount => hash that contains a counter for each feature type. It is used to create uniq ID
# 		example: $miscCount->{$primary_tag}++;
# $uniqID => hash of uniqID (UniqID link to the original ID)
# 		example: $uniqID->{$uID}=$id;
# $uniqIDtoType => hash to keep track about the feature type linked to the uniqID (useful for handling SPREADFEATURES)
# 		example: $uniqIDtoType->{$uID}=$primary_tag;
# $locusTAG_uniq => hash of comon tag found when reading the grouped features sequentialy
# 		example: $locusTAG_uniq->{'level1'}{$id}=$id;
# $infoSequential => hash that contains features grouped together in a sequential order. (Useful as example when no Parent tag and/or locus tag missing)
# 		structure at level1: $infoSequential->{$id}{'level1'}=$id;
# 		structure at other levels: $infoSequential->{$locusTAGvalue}{lc($l2_id)}{'level3'}}, [$feature1,$feature2] ;
# $comonTagAttribute => List of tags that could be use to gather features together (related to a same locus). Useful to keep track of features from a same locus that are spread whithin the file
# 		example: [locus_tag, gene_id]
# $last_locusTAGvalue => String: Last locus tag that has been met when parsing the file (If no locus tag found it will be the last feature L1 ID)
# 		example: CLUHARG00000008717
# $last_l1_f => String: Last l1 feature that has been met when parsing the file
# 		example: scaffold625	maker	gene	341518	341628	.	+	.	ID=CLUHARG00000008717
# $last_l2_f => String: Last L2 feature that has been met when parsing the file
# 		example: scaffold625	maker	mRNA	341518	341628	.	+	.	ID=CLUHART00000008717;Parent=CLUHARG00000008717
# $last_l3_f => String: Last L3 feature that has been met when parsing the file
# 		example: scaffold625	maker	exon	341518	341628	.	+	.	ID=CLUHART00000008717:exon:1406;Parent=CLUHART00000008717
# $verbose => INT: Verbose level. Bigger the value is, deeper the information sould be.
# 		example: 2
# ====== OUTPUT======= : Omniscient Hash
sub manage_one_feature{
	
	my ($ontology, $feature, $omniscient, $mRNAGeneLink, $duplicate, $miscCount, $uniqID, $uniqIDtoType, $locusTAG_uniq, $infoSequential, $comonTagAttribute, $last_locusTAGvalue, $last_l1_f, $last_l2_f, $last_l3_f, $last_f, $lastL1_new, $verbose)=@_;

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

#		+-------------------------------+
#		|	CKECK SEQUENCE ONTOLOGY		|
#		+-------------------------------+
		if(! exists($ontology->{$primary_tag}) ) {
			warn "GLOBAL@"."ontology1@".$primary_tag."@";
		}

#	+--------------------------------------------------------+
#		+-----------------------------------------------+
#		|	MANAGE LEVEL1 => feature WITHOUT parent		|
#		+-----------------------------------------------+
#	+--------------------------------------------------------+
	    if( get_level($feature) eq 'level1' ) {

	    	##########
			# get ID #
			#my $had_ID=undef;
	    	#if($feature->has_tag('ID')){$had_ID;}
	    	$id = lc(_check_uniq_id($omniscient, $miscCount, $uniqID, $uniqIDtoType, $feature));

	    	#####################
	    	# Ckeck duplication #
    		if(! _it_is_duplication($duplicate, $omniscient, $uniqID, $feature)){
    			
	    		################
	    		# Save feature #
	    		$last_l1_f = $feature;
	    		print "Push-L1-omniscient level1 || $primary_tag || $id = ".$feature->gff_string()."\n" if ($verbose >=2);
    			$omniscient->{"level1"}{$primary_tag}{$id}=$feature;
 	       		$locusTAG_uniq->{'level1'}{$id}=$id;

 	       		#############
 	       		# COMON TAG #
		        $locusTAGvalue =_get_comon_tag_value(\@comonTagList, $feature, $locusTAG_uniq, 'level1');

				if($locusTAGvalue){
					print "Push-L1-sequential $locusTAGvalue || level1 == $id\n" if ($verbose >=2);
					$locusTAG_uniq->{'level1'}{$locusTAGvalue}=$id;
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

		    return $id, $last_l1_f, $last_l2_f, $last_l3_f, $last_l1_f, $lastL1_new;
	    }

#	+-----------------------------------------------------------------------+
#		+---------------------------------------------------------------+
#		|	MANAGE LEVEL2 => feature WITH child and WITH parent		    |
#		+---------------------------------------------------------------+
#	+-----------------------------------------------------------------------+
      	elsif ( get_level($feature) eq 'level2' ) {

    		#reinitialization
    		$last_l3_f=undef;

    		##########
			# get ID #
	    	$id = _check_uniq_id($omniscient, $miscCount, $uniqID, $uniqIDtoType, $feature);

			##############
			# get Parent #
			if($feature->has_tag('Parent')){
				$parent = lc($feature->_tag_value('Parent'));
				$locusTAGvalue=$parent;
			}
			elsif($feature->has_tag('gene_id') ){
				$parent = lc($feature->_tag_value('gene_id'));
				create_or_replace_tag($feature,'Parent',$feature->_tag_value('gene_id')); #modify Parent To keep only one
				$locusTAGvalue=$parent;
			}
			else{warn "WARNING gff3 reader level2 : No Parent attribute found for @ the feature: ".$feature->gff_string()."\n";

				#################
 	       		# COMON TAG PART1
				$locusTAGvalue =_get_comon_tag_value(\@comonTagList, $feature, $locusTAG_uniq, 'level2');


				######################
				# NEED THE LEVEL1 ID #
				my $l1_ID="";
				# If I don't have a last_l1_f I create one. The Id can be used as comonTag. 
				# The feature can also be used later if comon_tag was existing, but mising for one of the feature. 
				# If we have comon tag, we check we are changing from the previous one before to create a new level1 feature. 
				# It's to deal with potential level2 (like mRNA isoforms).
				if(! $last_l1_f or ($locusTAGvalue and ($locusTAGvalue ne $last_locusTAGvalue) ) ){
					print "create L1 feature\n" if ($verbose >=3);
					$l1_ID = _create_ID($miscCount, $uniqID, $uniqIDtoType, $primary_tag, $id, "nbis_noL1id");					
					$last_l1_f = clone($feature);
					create_or_replace_tag($last_l1_f,'ID',$l1_ID); #modify Parent To keep only one
					$last_l1_f->primary_tag('gene');
					$lastL1_new = 1;
				}
				else{ # case where previous level1 exists
					# Stricly sequential at level2 feature. We create a new L1 at every L2 met except if two L2 are in a row
					if ( ($lastL1_new and not exists(LEVEL2->{$last_f->primary_tag}) ) and (!$locusTAGvalue or ($locusTAGvalue ne $last_locusTAGvalue) ) ){ # if previous L1 newly created and last feature is not f2 (if several f2 in a row we attach them to the same newly created l1 feature)
						print "create L1 feature stritcly\n" if ($verbose >=3);
						$l1_ID = _create_ID($miscCount, $uniqID, $uniqIDtoType, $primary_tag, $id, "nbis_noL1id");					
						$last_l1_f = clone($feature);
						create_or_replace_tag($last_l1_f,'ID',$l1_ID); #modify Parent To keep only one
						$last_l1_f->primary_tag('gene');
						$lastL1_new = 1;
						$lastL1_new = 1;
					}
					else{
						print "take last L1 feature\n" if ($verbose >=3); 
						$l1_ID=$last_l1_f->_tag_value('ID');
						$lastL1_new = undef;
					}
				}
				create_or_replace_tag($feature,'Parent',$l1_ID); #modify Parent To keep only one

				#################
 	       		# COMON TAG PART2
				if($locusTAGvalue){ #Previous Level up feature had a comon tag
					print "Push-L2-Sequential-1 $locusTAGvalue || ".lc($id)." || level2 == ".$feature->gff_string."\n" if ($verbose >=2);
		    		$infoSequential->{$locusTAGvalue}{lc($id)}{'level2'} = $feature ;
		
			    	return $locusTAGvalue, $last_l1_f, $feature, $last_l3_f, $feature, $lastL1_new;								#### STOP HERE AND RETURN
				}
				else{					
					
					print "Push-L2-omniscient-2: level2 || ".$primary_tag." || ".lc($l1_ID)." == ".$feature->gff_string."\n" if ($verbose >=2);
					push (@{$omniscient->{"level2"}{$primary_tag}{lc($l1_ID)}}, $feature);
					
					# keep track of link between level2->leve1 #
	  				if (! exists ($mRNAGeneLink->{lc($id)})){ 
						$mRNAGeneLink->{lc($id)}=$l1_ID;
		 			}	

					return lc($l1_ID) , $last_l1_f, $feature, $last_l3_f, $feature, $lastL1_new;								#### STOP HERE AND RETURN
				}
			}

			#####################
	    	# Ckeck duplication #
    		if(! _it_is_duplication($duplicate, $omniscient, $uniqID, $feature)){

				############################################
				# keep track of link between level2->leve1 #
	  			if (! exists ($mRNAGeneLink->{lc($id)})){ 
					$mRNAGeneLink->{lc($id)}=$parent;
		 		}
		 		
		 		####################
		 		# SAVE THE FEATURE #
		 		print "Push-L2-omniscient-3 level2 || $primary_tag || $parent == feature\n" if ($verbose >=2);
	      		push (@{$omniscient->{"level2"}{$primary_tag}{lc($parent)}}, $feature);
	      	}
	      	return $last_locusTAGvalue, $last_l1_f, $feature, $last_l3_f, $feature, $lastL1_new;
      	}

#	+--------------------------------------------------------+
#		+-----------------------------------------------+
#		|	MANAGE LEVEL3 => feature WITHOUT child 		|
#		+-----------------------------------------------+
#	+--------------------------------------------------------+
      	elsif ( get_level($feature) eq 'level3' ){
      		
      		##########
			# get ID #
	    	$id = _check_uniq_id($omniscient, $miscCount, $uniqID, $uniqIDtoType, $feature);
	    	
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
				$locusTAGvalue =_get_comon_tag_value(\@comonTagList, $feature, $locusTAG_uniq, 'level3');

				######################
				# NEED THE LEVEL2 ID #
				my $l2_id="";
				
				#To keep track of locus tag that has been spread over the file, and a piece is found later
				my $skip_last_l2=undef;
				if($last_l2_f and $locusTAGvalue){
					if(exists_keys ($locusTAG_uniq, ('linkl2l1', lc($last_l2_f->_tag_value('ID') ) ) ) ){
						if (lc($locusTAG_uniq->{'linkl2l1'}{lc( $last_l2_f->_tag_value('ID') )}) ne lc($locusTAGvalue)){
							$skip_last_l2=1;
						}
					}
				}

				# Just to avoid to have parent undef in case there is no parent feature define for the last_l2_f 
				my $parent_of_last_l2 = "@@@@";
				if($last_l2_f and $last_l2_f->has_tag('Parent')){ $parent_of_last_l2 = lc($last_l2_f->_tag_value('Parent')); }


				# case where No level2 feature defined yet - I will need a bucketL2 
				# OR comon tag changed (= level1/level2 different) so we have to create a new level2 tag - but only if the last_comon tag is different as the parent of the last_l2_f (In that case we can use the last L2 feature. It was missing the comon tag in it). 
				if(! $last_l2_f or ($locusTAGvalue and ($locusTAGvalue ne $last_locusTAGvalue) and $last_locusTAGvalue ne $parent_of_last_l2 or $skip_last_l2)  ){
					print "Create L2 feature $locusTAGvalue ne $last_locusTAGvalue !\n" if ($verbose >=3);
					
					#######################
					# Change referenrtiel => based on the last L2 link to this locus
					#######################
					if(exists_keys ($locusTAG_uniq, ('level2',$locusTAGvalue) ) ){ #first of

						$last_l2_f = @{$locusTAG_uniq->{'level2'}{$locusTAGvalue}}[$#{$locusTAG_uniq->{'level2'}{$locusTAGvalue}}]; # case were locus already met before (feature are spred within the file), we link the L3 to the last l2 of this locus.
						$l2_id = $last_l2_f->_tag_value('ID');
						foreach my $tag_l1 (keys %{$omniscient->{'level1'}}){
							if(exists_keys ($omniscient,('level1', $tag_l1, $locusTAG_uniq->{'linkl2l1'}{lc($l2_id)}))){
								$last_l1_f = $omniscient->{'level1'}{$tag_l1}{$locusTAG_uniq->{'linkl2l1'}{lc($l2_id)}};
							}
						}
					}
					else{

						$l2_id = _create_ID($miscCount, $uniqID, $uniqIDtoType, $primary_tag, $id, PREFIXL2);					
						$last_l2_f = clone($feature);
						create_or_replace_tag($last_l2_f,'ID',$l2_id); #modify Parent To keep only one
						$last_l2_f->primary_tag('RNA');
					}
				}
				else{ # case where previous level2 exists
					$l2_id=$last_l2_f->_tag_value('ID');
				}

				create_or_replace_tag($feature,'Parent',$l2_id); #modify Parent To keep only one

				#############
 	       		# COMON TAG  Part2
				if($locusTAGvalue){ #Previous Level up feature had a comon tag					
					print "Push-L3-sequential-1 $locusTAGvalue || ".lc($l2_id)." || level3 == ".$feature->gff_string."\n" if ($verbose >=2);
					### TAKE LAST L2 of the locus tag iF exist !
		    		push( @{$infoSequential->{$locusTAGvalue}{lc($l2_id)}{'level3'}}, $feature );
			    	return $locusTAGvalue, $last_l1_f, $last_l2_f, $feature, $feature, $lastL1_new;								#### STOP HERE AND RETURN
				}
				else{# No comon tag found 
					######################
					# NEED THE LEVEL1 ID #
					if(!$last_l1_f and $last_l3_f){ #particular case : Two l3 that follow each other, but first one has locus_tag but not the second
						print "Push-L3-sequential-2 $last_locusTAGvalue || ".lc($l2_id)." || level3 == ".$feature->gff_string."\n" if ($verbose >=2);
						push( @{$infoSequential->{$last_locusTAGvalue}{lc($l2_id)}{'level3'}}, $feature );
						return $last_locusTAGvalue, $last_l1_f, $last_l2_f, $feature, $feature, $lastL1_new;			
					}
					else{
						my $l1_id="";			
						if($last_l1_f){ # case where previous level1 exists
							$l1_id=$last_l1_f->_tag_value('ID');
						}
						else{ # case where No level1 feature defined yet - I will need a bucketL1
							$l1_id = _create_ID($miscCount, $uniqID, $uniqIDtoType, $primary_tag, $id, "nbis_noL1id");					
							$last_l1_f = clone($feature);
							create_or_replace_tag($last_l1_f,'ID',$l1_id); #modify Parent To keep only one
							$last_l1_f->primary_tag('gene');
						}

						#push( @{$infoSequential->{lc($l1_id)}{lc($l2_id)}{'level3'}}, $feature );
						print "Push-L3-omiscient-3: level3 ".$primary_tag." || ".lc($l2_id)." == ".$feature->gff_string."\n" if ($verbose >=2);
						push (@{$omniscient->{"level3"}{$primary_tag}{lc($l2_id)}}, $feature);
						return $l2_id, $last_l1_f, $last_l2_f, $feature, $feature, $lastL1_new;								#### STOP HERE AND RETURN
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
							_check_uniq_id($omniscient, $miscCount, $uniqID, $uniqIDtoType, $feature_clone); #Will change the ID if needed

							print "Push-L3-omniscient-4 level3 || $primary_tag || ".lc($parent)." == feature_clone\n" if ($verbose >=2);
							push (@{$omniscient->{"level3"}{$primary_tag}{lc($parent)}}, $feature_clone);
						}
						
						# It is the first parent. Do not clone the feature
						else{
							create_or_replace_tag($feature,'Parent',$parent); #modify Parent To keep only one
							print "Push-L3-omniscient-5 level3 || $primary_tag || ".lc($parent)." == feature\n" if ($verbose >=2);
							push (@{$omniscient->{"level3"}{$primary_tag}{lc($parent)}}, $feature);
						}
					}
					else{ #the simpliest case. One parent only
						print "Push-L3-omniscient-6 level3 || $primary_tag || ".lc($parent)." == feature\n".$feature->gff_string."\n" if ($verbose >=2);
						push (@{$omniscient->{"level3"}{$primary_tag}{lc($parent)}}, $feature);
					}
				}

				#Level3 key exists
				else{
					
					# It is a multiple parent case # Not the first parent, we have to clone the feature !!
					if($cptParent > 1){ #several parent, and we are not looking the first one
						
						my $feature_clone=clone($feature);
						create_or_replace_tag($feature_clone,'Parent',$parent); #modify Parent To keep only one
						_check_uniq_id($omniscient, $miscCount, $uniqID, $uniqIDtoType, $feature_clone); #Will change the ID if needed

						if( ! _it_is_duplication($duplicate, $omniscient, $uniqID, $feature_clone) ){
							print "Push-L3-omniscient-8 level3 || $primary_tag || ".lc($parent)." == feature_clone\n" if ($verbose >=2);
							push (@{$omniscient->{"level3"}{$primary_tag}{lc($parent)}}, $feature_clone);
						}
								
					}
					elsif($allParent > 1){ # It is a multiple parent case #several parent, but we are looking the first one

						# It is the first parent. Do not clone the feature
						create_or_replace_tag($feature,'Parent',$parent); #modify Parent To keep only one	
						if( ! _it_is_duplication($duplicate, $omniscient, $uniqID, $feature) ){
							print "Push-L3-omniscient-9 level3 || $primary_tag || ".lc($parent)." == feature\n" if ($verbose >=2);
							push (@{$omniscient->{"level3"}{$primary_tag}{lc($parent)}}, $feature);
						}

					}
					#It is not a duplicated feature => save it in omniscient
					elsif( ! _it_is_duplication($duplicate, $omniscient, $uniqID, $feature) ){
 						#the simpliest case. One parent only
						print "Push-L3-omniscient-10 level3 || $primary_tag || ".lc($parent)." == feature\n".$feature->gff_string."\n" if ($verbose >=2);
						push (@{$omniscient->{"level3"}{$primary_tag}{lc($parent)}}, $feature);
					}
				}
      		}
      		return $last_locusTAGvalue, $last_l1_f, $last_l2_f, $feature, $feature, $lastL1_new;
      	}


#	+--------------------------------------------------------+
#		+-----------------------------------------------+
#		|	MANAGE THE REST => feature UNKNOWN   		| # FEATURE NIT DEFINE IN ANY OF THE 3 LEVELS
#		+-----------------------------------------------+
#	+--------------------------------------------------------+ 	
      	else{
      		warn "gff3 reader warning: primary_tag error @ ".$primary_tag."still not taken in account ! Please modify the code to define on of the three level of this feature.\n";
      		warn "GLOBAL@"."parser1@".$primary_tag."@";
      		return $locusTAGvalue, $last_l1_f, $last_l2_f, $last_l3_f, $feature, $lastL1_new;
      	}	

    print "Congratulation ! Read this line is not normal !! Please contact the developer !!!\n";exit;
}

##==============================================================

#check if the comom tag is present among the attributes
# return lower case value
sub _get_comon_tag_value{
	my ($comonTagList, $feature, $locusTAG_uniq, $level)=@_;

	my $locusName=undef;

   	foreach my $tag (@$comonTagList){

   		#check if we hace the tag
		if($feature->has_tag($tag)){
		    $locusName=lc($feature->_tag_value($tag)); #get the value
		    
		    if(exists_keys ($locusTAG_uniq, ('level1',$locusName) ) ){
		    	$locusName = $locusTAG_uniq->{'level1'}{$locusName};
		    	last;
		    }
		    else{
		    	$locusTAG_uniq->{'level1'}{$locusName} = $locusName; #save it
		    	last;
			}
		}
	}

	if($level eq 'level2' and $locusName){
		if(! exists_keys ($locusTAG_uniq, ('level2',$locusName, lc($feature->_tag_value('ID'))) ) ){
			push @{$locusTAG_uniq->{'level2'}{$locusName}}, $feature;
			$locusTAG_uniq->{'linkl2l1'}{lc($feature->_tag_value('ID'))} =  $locusName;
		}
	}

	# In case where no parent, no comon tag, and no sequential, we cannot deal at all with it !!!!
	if(! $locusName and $level ne 'level1'){
		warn "WARNING gff3 reader: Hmmm, be aware that your feature doesn't contain any Parent and locus tag. No worries, we will handle it by considered it as striclty sequential. If you are not ok with that, provide an ID or a comon tag by locus. @ the feature is:\n".$feature->gff_string()."\n";
	}

	return $locusName;
}

#feature is not yet saved in omniscient !
sub _it_is_duplication{
	my ($duplicate, $omniscient, $uniqID, $feature)=@_;

	my $is_dupli=undef;
	my $potentialList=undef;

	my $level = get_level($feature);
	my $primary_tag = lc($feature->primary_tag);

	my $id = $uniqID->{$feature->_tag_value('ID')}; # check the original ID

	if($level eq "level1"){
		if(! exists_keys($omniscient,($level, $primary_tag, lc($id) ))){
			return $is_dupli; #return is not a dupli
		}
		else{
			$potentialList=$omniscient->{$level}{$primary_tag}{lc($id)}; #push the feature L1 in potentialList
		}
	}
	else{ #feature l2 or l3

		my @parent = $feature->get_tag_values('Parent');
		foreach my $one_parent_ID (@parent){

			my $one_parent_uID = $one_parent_ID; # In case where the parent have not yet been processed, we cannot have his uID, we will check the current ID
			if ( exists_keys($uniqID, ($one_parent_ID)) ){
				$one_parent_uID = lc($uniqID->{$one_parent_ID}); # check the original ID
			}

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
sub get_level{
	my ($feature)=@_;

	my $source_tag = lc($feature->source_tag);	
	my $primary_tag = lc($feature->primary_tag);

	my $level=undef;

	#########################################
	## PECULIARITIES FROM HAVANA / ENSEMBL ##
	if ($source_tag eq "ensembl" ){
		if ( $primary_tag eq "rna" ) {return 'level1';} #particularity ENSEMBL
	}
	if ( ($source_tag =~ "havana" or $source_tag =~ "ensembl") and ($primary_tag eq  "processed_transcript" or $primary_tag eq  "pseudogene" ) ){ #By default processed_transcript is l2 and pseudogene is l1
		if ($feature->has_tag('Parent')){return "level2" ;}
		else{return "level1" ;}
	}
	## PECULIARITIES FROM HAVANA / ENSEMBL ##
	#########################################

	if (exists(LEVEL1->{$primary_tag}) ){
		return 'level1';
	}
	elsif(exists(LEVEL2->{$primary_tag}) ){
		return 'level2';
	}
	elsif(exists(LEVEL3->{$primary_tag}) ){
		return 'level3';
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
	my	($omniscient, $miscCount, $uniqID, $uniqIDtoType, $feature)=@_;
	
	my $uID=undef;
	my $primary_tag = lc($feature->primary_tag);

	my $id=undef;
	if($feature->has_tag('ID')){ #has the tag
		$id = $feature->_tag_value('ID');
	}
	elsif($feature->has_tag($primary_tag."_id") ){
		$id = $feature->_tag_value($primary_tag."_id");
		create_or_replace_tag($feature, 'ID', $id);
	}

	# CHECK THE ID TO SEE IF IT's uniq, otherwise we have to create a new uniq ID
	if($id){
		# In case of non-spreadfeature (avoid CDS and UTR that can share identical IDs)
		if(! SPREADFEATURE->{$primary_tag}){ 
			$uID = _create_ID($miscCount, $uniqID, $uniqIDtoType, $primary_tag, $id, 'nbis_NEW'); #method will push the uID
			if(	$id ne $uID ){ #push the new ID if there is one
				create_or_replace_tag($feature, 'ID', $uID);
			}
		}
		# In case of spreadfeature ( CDS and UTR that can share identical IDs)
		else{
			# First time we see this ID => No problem;
	 		if(! exists($uniqID->{$id})){	
			 	#push the uID
			 	$uID = $id;
			 	$uniqID->{$uID}=$id;
			 	$uniqIDtoType->{$id}=$primary_tag;
			}
		# NOT the first time we have this ID	
			# check if it's the same type (To not mix a same ID between UTR and CDS);
			elsif( $uniqIDtoType->{$id} eq $primary_tag ){ # Same type, so we can keep this ID, let's continue
			 	$uID = $id;
			}
			else{ # The spreadfeature type is different
				# Let's check if one of the same type is already in omniscient (THE ID could be linked to a non-spreadfeature), in that case we keep the ID already given.
				if( $feature->has_tag('Parent') ){
					if ( exists_keys( $omniscient, ('level3', $primary_tag, lc($feature->_tag_value('Parent')) ) ) ){
						$uID = 	@{ $omniscient->{'level3'}{$primary_tag}{lc($feature->_tag_value('Parent'))} }[0]->_tag_value('ID');
					}
				}
				if(! $uID){ #ID already taken by another feature type, and we do not have ID already existing of this feature type within omniscient, let's create a new ID
					$uID = _create_ID($miscCount, $uniqID, $uniqIDtoType, $primary_tag, $id, 'nbis_NEW'); #method will push the uID 	
				}
				if(	$id ne $uID ){ #push the new ID if there is one
				 		create_or_replace_tag($feature, 'ID', $uID);
				}
			}
		}
	}
	else{ #tag absent
		my $level = get_level($feature);
		if($level ne 'level3'){
			warn "gff3 reader error ".$level .": No ID attribute found @ for the feature: ".$feature->gff_string()."\n";
		}
		$miscCount->{$primary_tag}++;
		$id = $primary_tag."-".$miscCount->{$primary_tag}; # create an ID and then check if not already taken
		$uID = _create_ID($miscCount, $uniqID, $uniqIDtoType, $primary_tag, $id, 'nbis_NEW'); #method will push the uID
		create_or_replace_tag($feature, 'ID', $uID);
	}

	return $uID;
}

# create the ID and add it to the feature.
sub _create_ID{
	my	($miscCount, $uniqID, $uniqIDtoType, $primary_tag, $id, $prefix)=@_;
	
	my $key;

	if($prefix){
		$key=$prefix."-".$primary_tag;
	}
	else{
		$key=$primary_tag;
	}

	my $uID=$id;
	while( exists_keys($uniqID, ($uID) )){	 #loop until we found an uniq tag	
		$miscCount->{$key}++;
		$uID = $key."-".$miscCount->{$key};
	}

	#push the new ID	
	$uniqID->{$uID}=$id;
	$uniqIDtoType->{$uID}=$primary_tag;

	return $uID;
}

# check if mRNA have is Parental gene existing. If not we create it.
sub _check_gene_link_to_mrna{
	my ($hash_omniscient, $verbose)=@_;
	my $resume_case=undef;

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
				$resume_case++;
				print "WARNING gff3 reader level2 : No Parent feature found with the ID @ ".$id_l1.". We will create one.\n" if ($verbose >= 2);
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
			}
		}
	}

	print "We create $resume_case level1 features that were missing\n" if($verbose >= 1 and $resume_case);
}

# @Purpose: Remove the level1 feature that havn't subfeature linked to it. Before to remove it check if L3 is linked to it. In that case it is a format error that we will fix !
# @input: 1 => hash(omniscient hash)
# @output: none
sub _remove_orphan_l1{
	my ($hash_omniscient, $miscCount, $uniqID, $uniqIDtoType, $mRNAGeneLink, $verbose)=@_;
	my $resume_case=undef;

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
     			if ( ! refseq_case($hash_omniscient, $hash_omniscient->{'level1'}{$tag_l1}{$id_l1}, $miscCount, $uniqID, $uniqIDtoType, $mRNAGeneLink, $verbose)){
     				warn "WARNING gff3 reader level1 : No child feature found for the $tag_l1 with ID @ ".$id_l1.".\n";
     				print "One orphan removed \n" if ($verbose >= 2 );
     				$resume_case++;
     			}
 			    delete $hash_omniscient->{'level1'}{$tag_l1}{$id_l1}; # delete level1 // In case of refseq the thin has been cloned and modified, it is why we nevertheless remove it
		    }
 	 	}
 	}
 	print "We create $resume_case level1 features that were missing\n" if($verbose >= 1 and $resume_case);
}

# Level3 related to level1 wihtout level2 defined
# so the id of L1 is transferred to l2 to keep the parent attribute of l3 correct
sub refseq_case{ # refseq case as example
 	my ($omniscient, $l1_feature, $miscCount, $uniqID ,$uniqIDtoType, $mRNAGeneLink, $verbose)=@_;

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
 			        $l1_ID = _create_ID($miscCount, $uniqID, $uniqIDtoType, lc($l1_feature->primary_tag), $l1_ID, 'nbis_NEW');
 				  	create_or_replace_tag($l1_feature,'ID', $l1_ID); #modify ID to replace by parent value
 				  	my $l1_feature_clone = clone($l1_feature);#create a copy of the first mRNA feature;
 				  	$omniscient->{"level1"}{lc($l1_feature->primary_tag)}{lc($l1_ID)} = $l1_feature_clone;

 				  	#Modify parent L2
 				  	create_or_replace_tag($l2_feature,'Parent', $l1_ID); #modify ID to replace by parent value
 				  	push(@{$omniscient->{"level2"}{lc($l2_feature->primary_tag)}{lc($l1_ID)}}, $l2_feature);

 				  	#fill the $mRNAGeneLink hash
 				  	$mRNAGeneLink->{lc($id_l2)} = $l1_ID; # Always need to keep track about l2->l1, else the method _check_l3_link_to_l2 will recreate a l1 thinking this relationship is not fill
 				  	print "L1 child was directly linked to L3. Corrected by removing the L1 and creating new L1 and L2 features\n" if($verbose >= 1);

 				  	return 1;
 		    }
 	 	}
 	}
 	return 0;
}

# @Purpose: Check relationship betwwen L3 and L2. If L2 is missing we create it. When creating L2 missing we create as well L1 if missing too.
# @input: 4 => hash(omniscient hash), hash(mRNAGeneLink hash), hash(miscCount hash), hash(uniqID hash)
# @output: none
sub _check_l3_link_to_l2{
	my ($hash_omniscient, $mRNAGeneLink, $miscCount, $uniqID, $uniqIDtoType, $verbose)=@_;
	my $resume_case=undef;

 	foreach my $tag_l3 (keys %{$hash_omniscient->{'level3'}}){

 	  	foreach my $id_l2 (keys %{$hash_omniscient->{'level3'}{$tag_l3}}){

 	  		#check if L2 exits
 	  		if (! exists($mRNAGeneLink->{lc($id_l2)}) ) {
 	  			$resume_case++;
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

				my $new_ID_l1 = _check_uniq_id($hash_omniscient, $miscCount, $uniqID, $uniqIDtoType, $l1_feature);
				create_or_replace_tag($l1_feature,'ID', $new_ID_l1); #modify ID to replace by parent value
				
			#finish fill Level2
				create_or_replace_tag($l2_feature, 'Parent', $new_ID_l1); # remove parent ID because, none.
				#save new feature L2
				push (@{$hash_omniscient->{"level2"}{$primary_tag_l2}{lc($new_ID_l1)}}, $l2_feature);

			#finish fill Level1
				check_level1_positions($hash_omniscient, $l1_feature);	# check start stop if isoforms exists
				#save new feature L1
				$hash_omniscient->{"level1"}{$primary_tag_l1}{lc($new_ID_l1)} = $l1_feature; # now save it in omniscient
				$mRNAGeneLink->{lc($id_l2)} = $new_ID_l1;

				print "L1 and L2 created, \n" if($verbose >= 1);
 	  		}
 	  	}
 	}
 	print "We fixed $resume_case cases where L2 and L1 features were missing\n" if($verbose >= 1 and $resume_case);	 
}

# @Purpose: Check L3 features. If exon are missing we create them. We go through all features of level3 and check them by type, if two should be merged, we do it (CDS 1-50 and 51-100, must be CDS 1-100).
# @input: 3 =>  hash(omniscient hash), hash(miscCount hash), hash(uniqID hash)
# @output: none
sub _check_exons{
	my ($hash_omniscient, $mRNAGeneLink, $miscCount, $uniqID, $uniqIDtoType, $verbose)=@_;
	my $resume_case=undef; my $resume_case2=undef; my $resume_case3=undef;

	my %checked;
	foreach my $tag_l3 (keys %{$hash_omniscient->{'level3'}}){
		if ($tag_l3 ne "exon"){
 	  		foreach my $id_l2 (keys %{$hash_omniscient->{'level3'}{$tag_l3}}){
 	  			
 	  			if( ! exists_keys(\%checked,($id_l2)) ){ #l2 already checked

 	  				my $feature_example=undef; # will be used to create the exon features
	 	  			my $list_location_Exon=[];
	 	  			my $list_location_NoExon=[];
	 	  			my $list_location_NoExon_overlap=[];
	 	  			my %createIT; # will be usefull to list the feature to create
#				 	+-----------------------------------------------------
#					| 			Go through l3 and save info needed		 |
#				 	+-----------------------------------------------------

	 	  			foreach my $tag_l3 (keys %{$hash_omniscient->{'level3'}}){
	 	  				
	 	  				# LIST NON-EXON LOCATIONS THAT NEED TO BE IN AN EXON LOCATION
	 	  				if ($tag_l3 ne "exon" and LEVEL3->{$tag_l3} eq "exon" ){

				 	  		if( exists_keys($hash_omniscient,('level3',$tag_l3, $id_l2)) ){				 	  		

				 	  			my $list_location_l3=[];
				 	  			foreach my $l3_feature (@{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2}}){

				 	  				if(! $feature_example){
				 	  					$feature_example=$l3_feature;
				 	  				}
				 	  				
				 	  				my $locationRefList=[[[$l3_feature->_tag_value('ID')] ,int($l3_feature->start), int($l3_feature->end)]];
				 	  				$list_location_l3 = _manage_location($locationRefList, $list_location_l3, 'adjacent', 0); # we check first in overlap mode to check if badly define features exists 				
				 	  			}

				 	  			#Rare case when a feature of the same type is badly defined
				 	  			if($#{$list_location_l3} < $#{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2}}){
				 	  				warn "Peculiar rare case, we should remove existing $tag_l3 which are supernumerary.Not yet implemented\n";
				 	  			}
				 	  			push @{$list_location_NoExon_overlap}, @{$list_location_l3}; #list of all feature that has been checked in overlap mode
				 	  		}
				 	  	}

				 	  	# LIST EXON LOCATIONS
				 	  	elsif($tag_l3 eq "exon"){

							if( exists_keys($hash_omniscient,('level3',$tag_l3, $id_l2)) ){				 	  		

				 	  			foreach my $l3_feature (@{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2}}){

				 	  				if(! $feature_example){
				 	  					$feature_example=$l3_feature;
				 	  				}
				 	  				my $locationRefList=[[[$l3_feature->_tag_value('ID')] ,int($l3_feature->start), int($l3_feature->end)]];
				 	  				$list_location_Exon = _manage_location($locationRefList, $list_location_Exon , 'adjacent', 0);
				 	  			}
				 	  			
				 	  			#Rare case when a features are badly defined 
				 	  			# This approch works for exon because they have uniq ID
				 	  			if($#{$list_location_Exon} < $#{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2}}){
				 	  				warn "Peculiar rare case, we have to remove existing exon which are supernumerary. Parent is $id_l2\n";
				 	  				#remove the non needed features (they where wrong, and are unecessary)
									my @id_list2=();
									foreach my $locations (@{$list_location_Exon}){
										# If several value in the ID list, we have to avoid the first value (which is the one to keep), and remove all the other ones which are the Id to remove
										if(@{$locations->[0]} > 1){
											my $correct_ID = shift @{$locations->[0]};
											push @id_list2, @{$locations->[0]};
											@{$locations->[0]} = $correct_ID;
											foreach my $l3_feature (@{$hash_omniscient->{'level3'}{'exon'}{$id_l2} } ){
			 	  								if($l3_feature->_tag_value('ID') eq  $correct_ID){
			 	  									$l3_feature->start($locations->[1]);
			 	  									$l3_feature->end($locations->[2]);
			 	  								}
			 	  							}
										}
									}
				 	  				my @tag_list = ('all');
				 	  				my @id_list=($id_l2);
				 	  				$resume_case3 += $#id_list2;
				 	  				print "We remove the supernumerary @id_list2 exon(s)\n" if($verbose >= 2);
									remove_element_from_omniscient(\@id_list, \@id_list2, $hash_omniscient, 'level3', 'false', \@tag_list);
				 	  			}
				 	  		
				 	  		}				 	  		
				 	  	}
				 	}

				 	## check all NOn-exon in adjacent mater to have the correct list (allows to merge UTR and CDS to recreate exon )
				 	foreach my $location (@{$list_location_NoExon_overlap}){
				 		$list_location_NoExon = _manage_location([$location], $list_location_NoExon, 'adjacent', 0);
				 	}

				 	print "list_location_Exon: ".Dumper($list_location_Exon) if ($verbose >= 3); 
				 	print "list_location_NOEXON: ".Dumper($list_location_NoExon) if ($verbose >= 3);  
 	  				

#				 	+--------------------------------------------------------------------------------------------------------+
#					| 				COMPARE EXONS POSITION TO THOSE DESCRIBED BY NON-EXON FEATURES 						 |
#				 	+--------------------------------------------------------------------------------------------------------+

 	  				#Case where exon feature exists, we have to check them
	 	  			if( exists_keys($hash_omniscient,('level3','exon', $id_l2)) ){ #l3 feature but no exon among them... need to recreate them.

	 	  				#create string to comapre the 2 lists.
	 	  				my $list_location_Exon_joined="";
	 	  				foreach my $location (sort {$a->[1] <=> $b->[1] } @{$list_location_Exon}){
	 	  					$list_location_Exon_joined .= $location->[1].$location->[2];
	 	  				}
	 	  				my $list_location_NoExon_joined="";
	 	  				foreach my $location (sort {$a->[1] <=> $b->[1] } @{$list_location_NoExon}){
	 	  					$list_location_NoExon_joined .= $location->[1].$location->[2];
	 	  				}
	 	  				#If two lists are different we have to check/fix the difference
	 	  				# If no overlap we create the exon:
	 	  				# If overlap:  Redefine internal exon ; Redfine external exon only if too short.
						if($list_location_Exon_joined ne $list_location_NoExon_joined ){
							print "_check_exons EXON MISSING ! Let's check that !! \n" if ($verbose >= 2);

							my $location_cpt=0;
		 	  				foreach my $location (sort {$a->[1] <=> $b->[1] } @{$list_location_NoExon}){
		 	  					$location_cpt++;

		 	  					my $create_exon=1;
		 	  					my $new_location;
		 	  					my $overlap;
		 	  					
		 	  					foreach my $exon_location (sort {$a->[1] <=> $b->[1] } @{$list_location_Exon}){
		 	  						
		 	  						($new_location, $overlap) = _manage_location_lowLevel_adjacent($location, $exon_location); #there is an overlap if $new_location != $exon_location. If it's the same, we should check $overlap to be sure

		 	  						if($new_location->[1] < $exon_location->[1] or $new_location->[2] > $exon_location->[2] ){ #The exon_location has been modified by location... We have to remodelate the exon (only if fit some conditions) location to take the modification into account
			 	  						$create_exon=undef; # We must avoid to create exon because there is an overlap. 

		 	  							my $redefine_left=undef;
		 	  							my $redefine_right=undef;
		 	  							#first location => check left
			 	  						if($location_cpt == 1){
			 	  							if($new_location->[1] <  $exon_location->[1]){ $redefine_left = $new_location->[1];} # Modify only if it's more left
			 	  						}	
			 	  						#=> check left and right 
			 	  						if($location_cpt != 1 and $location_cpt != @$location){
			 	  							if($new_location->[1] <  $exon_location->[1]){ $redefine_left = $new_location->[1];}  # Modify only if it's more left
			 	  							if($new_location->[2] >  $exon_location->[2]){ $redefine_right = $new_location->[2];} # Modify only if it's more right
			 	  						}
			 	  						#last location => check right
			 	  						if($location_cpt == @$location){
			 	  							if($new_location->[2] >  $exon_location->[2]){ $redefine_right = $new_location->[2];} # Modify only if it's more right
			 	  						}

			 	  						foreach my $l3_feature (@{$hash_omniscient->{'level3'}{'exon'}{$id_l2} } ){
			 	  							if($l3_feature->_tag_value('ID') eq $exon_location->[0][0]){
			 	  								
			 	  								if($redefine_left){
			 	  									$l3_feature->start($new_location->[1]);
			 	  								}else{$redefine_left = $exon_location->[1];}
			 	  								
			 	  								if($redefine_right){
			 	  									$l3_feature->end($new_location->[2]);
			 	  								}else{$redefine_right = $exon_location->[2];}
			 	  								if($redefine_left or $redefine_right){$resume_case2++;}
				 	  							print "We modify the location of the existing exon !! ".$exon_location->[0][0]." ".$exon_location->[1]." ".$exon_location->[2]." to ".$redefine_left." ".$redefine_right."\n" if ($verbose >= 2);
			 	  								last;
			 	  							}
			 	  						}
			 	  					}
			 	  					elsif($overlap){ #location not modified but no overlap, so it means the exon is not defined !
			 	  						$create_exon=undef;
			 	  					}
			 	  				}

			 	  				if($create_exon){
			 	  					push @{$createIT{'exon'}}, $location;
			 	  				}
	 						}
	 					}
	 	  			}
	 	  			else{ $createIT{'exon'}=$list_location_NoExon;} # no exon exists, we have to create all of them
 
					# NOW CREATE EXON IF NECESSARY
					if(keys %createIT){
						foreach my $tag (keys %createIT){
					 	  	foreach my $location (@{$createIT{$tag}}){
					 	  		$resume_case++;
					 	  		my $feature_exon = clone($feature_example);#create a copy of a random feature l3;
								$feature_exon->start($location->[1]);
					 	  		$feature_exon->end($location->[2]);
					 	  		$feature_exon->frame(".");
					 	  		$feature_exon->primary_tag($tag);
					 	  		my $uID = _check_uniq_id($hash_omniscient, $miscCount, $uniqID, $uniqIDtoType, $feature_exon);
					 	  		create_or_replace_tag($feature_exon, 'ID', $uID); # remove parent ID because, none.
								#save new feature L2
								print "_check_exons Create one Exon : ".$feature_exon->gff_string."\n" if ($verbose >= 2);
								push (@{$hash_omniscient->{"level3"}{$tag}{$id_l2}}, $feature_exon);
					 	  	}
					 	}
				 	}

				 	#Check extremities of exons (If shorter we adapt to the mRNA size, else we adapt the L2 to the exon size)
	 	  			my $id_l1 = lc($mRNAGeneLink->{lc($id_l2)});
	 	  			my $getout=undef;
	 	  			foreach my $tag_l2 ( %{$hash_omniscient->{'level2'}} ){
	 	  				if( exists_keys($hash_omniscient,('level2', $tag_l2, $id_l1)) ){
	 	  					foreach my $l2_feature ( @{$hash_omniscient->{'level2'}{$tag_l2}{$id_l1}} ){
	 	  						if( lc($l2_feature->_tag_value('ID')) eq $id_l2 ){
	 	  						 	
	 	  						 	my $myLeftExtremity=$l2_feature->start();
	 	  						 	my $myRightExtremity=$l2_feature->end();
			 	  		
				 	  			 	my @list_exon = sort {$a->start <=> $b->start} @{$hash_omniscient->{'level3'}{'exon'}{$id_l2}};

				 	  			 	if( int($list_exon[0]->start) >  int($myLeftExtremity) ){
				 	  			 		print "_check_exons We modified the exon LEFT extremity from $id_l2! ".$list_exon[0]->start." <to> ".$myLeftExtremity."\n" if($verbose >= 1);;
				 	  			 		$list_exon[0]->start($myLeftExtremity);
				 	  			 	}
				 	  			 	if($list_exon[0]->start <  $myLeftExtremity){  #modify L2 
				 	  			 		$l2_feature->start($list_exon[0]->start);
				 	  			 		print "_check_exons We modified the L2 LEFT extremity !\n" if($verbose >= 1);
				 	  			 	}

				 	  			 	if($list_exon[$#list_exon]->end <  $myRightExtremity){
				 	  			 		print "_check_exons We modified the exon RIGHT extremity from $id_l2!".$list_exon[$#list_exon]->end." to ".$myRightExtremity."\n" if($verbose >= 1);
				 	  			 		$list_exon[$#list_exon]->end($myRightExtremity);  			 		
				 	  			 	}
				 	  			 	elsif($list_exon[$#list_exon]->end >  $myRightExtremity){ #modify L2 
				 	  			 		$l2_feature->end($list_exon[$#list_exon]->end);
				 	  			 		print "_check_exons We modified the L2 RIGHT extremity !\n" if($verbose >= 1);
				 	  			 	}

				 	  			 	$getout=1;
				 	  			 	last;
				 	  			}
	 	  					}
	 	  					if($getout){
	 	  						last;
	 	  					}
	 	  				}
	 	  			}

				 	#keep track of l2 checked (as we loop over L3, we meet several time the same l2)
		 	  		$checked{$id_l2}++;
	 	  		}
 	  		}
 	  	}
 	}
 	print "We create $resume_case exons that were missing\n" if($verbose >= 1 and $resume_case);
 	print "We modified $resume_case2 exons positions that were wrong\n" if($verbose >= 1 and $resume_case2);
 	print "We have removed $resume_case3 existing exons which were supernumerary\n" if($verbose >= 1 and $resume_case3);
}

# @Purpose: Check L3 features. If UTRS are missing we create them.
# @input: 3 =>  hash(omniscient hash), hash(miscCount hash), hash(uniqID hash)
# @output: none
sub _check_utrs{
	my ($hash_omniscient, $mRNAGeneLink, $miscCount, $uniqID, $uniqIDtoType, $verbose)=@_;
	my $resume_case=undef;my $resume_case2=undef;

	my %checked;
	foreach my $tag_l3 (keys %{$hash_omniscient->{'level3'}}){
		if ($tag_l3 ne "exon"){
 	  		foreach my $id_l2 (keys %{$hash_omniscient->{'level3'}{$tag_l3}}){
 	  			
 	  			if( ! exists_keys(\%checked,($id_l2)) ){ #l2 already checked

 	  				my $feature_example=undef; # will be used to create the exon features
	 	  			my $list_location_Exon=undef;
	 	  			my $list_location_CDS=undef;
	 	  			my $list_location_UTR=undef;
	 	  			
#				 	+-----------------------------------------------------
#					| 			Go through l3 and save info needed		 |
#				 	+-----------------------------------------------------

	 	  			foreach my $tag_l3 (keys %{$hash_omniscient->{'level3'}}){
	 	  				
				 	  	# LIST CDS LOCATIONS
	 	  				if ($tag_l3 eq "cds"){
				 	  		if( exists_keys($hash_omniscient,('level3',$tag_l3, $id_l2)) ){				 	  		

				 	  			foreach my $l3_feature (@{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2}}){

				 	  				my $locationRefList=[[[$l3_feature->_tag_value('ID')] ,int($l3_feature->start), int($l3_feature->end)]];
				 	  				$list_location_CDS = _manage_location($locationRefList, $list_location_CDS, 'adjacent', $verbose);
				 	  			}
				 	  		}
				 	  	}

				 	  	# LIST UTR LOCATIONS
	 	  				if ($tag_l3 =~ "utr"){
				 	  		if( exists_keys($hash_omniscient,('level3',$tag_l3, $id_l2)) ){				 	  		

				 	  			foreach my $l3_feature (@{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2}}){

				 	  				my $locationRefList=[[[$l3_feature->_tag_value('ID')] ,int($l3_feature->start), int($l3_feature->end)]];
				 	  				$list_location_UTR = _manage_location($locationRefList, $list_location_UTR, 'adjacent', $verbose);
				 	  			}
				 	  		}
				 	  	}

				 	  	# LIST EXON LOCATIONS
				 	  	elsif($tag_l3 eq "exon"){
							if( exists_keys($hash_omniscient,('level3',$tag_l3, $id_l2)) ){				 	  		

				 	  			foreach my $l3_feature (@{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2}}){

				 	  				if(! $feature_example){
				 	  					$feature_example=$l3_feature;
				 	  				}
				 	  				#print "exonFeature= ".$l3_feature->gff_string."\n";
				 	  				push @{$list_location_Exon}, [ [$l3_feature->_tag_value('ID')], int($l3_feature->start), int($l3_feature->end)] ;
				 	  			}
				 	  		}				 	  		
				 	  	}
				 	}

#				 	+-----------------------------------------------------
#					| 				HANDLE UTRs 						 |
#				 	+-----------------------------------------------------

					if( exists_keys($hash_omniscient,('level3','cds', $id_l2)) ){ #Check UTR only if CDS exists

						# Create list of UTR expected:
						my $list_location_UTR_expected=undef;
						my $expected_utr=1;
						
						foreach my $exon_location (sort {$a->[1] <=> $b->[1] } @{$list_location_Exon}){

			 	  			my $new_location;
			 	  			my $overlap;
			 	  			my $never_overlap="yes";
							foreach my $location_cds (sort {$a->[1] <=> $b->[1] } @{$list_location_CDS}){

								if( $location_cds->[1] > $exon_location->[2]){last;}
								if( $location_cds->[2] < $exon_location->[1]){next;}

								($new_location, $overlap) = _manage_location_lowLevel_inversed($location_cds, $exon_location, $verbose);

								if($overlap eq "perfect"){ $never_overlap=undef; $expected_utr=undef;last;}
								
								if($new_location->[1] != $exon_location->[1] and $new_location->[2] != $exon_location->[2] ){ #two UTR expected        =========================  exon
									print "creation utr push1\n" if($verbose >= 3);
									push @{$list_location_UTR_expected}, [undef, $exon_location->[1], $location_cds->[1]-1];				#								=======			CDS
									push @{$list_location_UTR_expected}, [undef, $location_cds->[2]+1, $exon_location->[2]];
									$never_overlap=undef;
									last;
								}	
								elsif($new_location->[1] != $exon_location->[1] or $new_location->[2] != $exon_location->[2] ){ #two UTR expected  {
									print "creation utr push2\n".Dumper($new_location)."\n" if($verbose >= 3);
									push @{$list_location_UTR_expected}, $new_location;
									$never_overlap=undef;
									last;
								}
							}
							if($never_overlap){ #case of UTR that match++ fully the exon
								push @{$list_location_UTR_expected}, $exon_location;
							}
						}

						print "list_location_UTR_expected: ".Dumper($list_location_UTR_expected) if ($verbose >= 3);
						print "list_location_UTR: ".Dumper($list_location_UTR) if ($verbose >= 3);  
 	  				
		 	  			# Compare UTR Present and UTR expected
	 	  				my $list_utr_to_create=undef;

	 	  				if($list_location_UTR){ #List UTR not empty
							if($list_location_UTR_expected){ #List UTR not empty
				 	  			foreach my $UTRexp_location (sort {$a->[1] <=> $b->[1] } @{$list_location_UTR_expected} ){
				 	  					
			 	  					my $create_utr=1;
			 	  					my $new_location;
			 	  					my $overlap;
			 	  					foreach my $UTR_location (sort {$a->[1] <=> $b->[1] } @{$list_location_UTR}){
			 	  						
			 	  						($new_location, $overlap) = _manage_location_lowLevel_inversed($UTR_location, $UTRexp_location, $verbose); #just to check that it overlaps

			 	  						if($overlap and ( $UTR_location->[1] != $UTRexp_location->[1] or $UTR_location->[2] != $UTRexp_location->[2] ) ){ #It overlaps and at least one location is different. We have to re-modelate the utr location to take the modification into account
				 	  						print "We modify the location of the existing utr: ".$UTR_location->[0][0]."  ".$UTR_location->[1]." ".$UTR_location->[2]." to ".$UTRexp_location->[1]." ".$UTRexp_location->[2]."\n" if ($verbose >= 3);
				 	  						$resume_case2++;
				 	  						$create_utr=undef;

				 	  						foreach my $tag_l3 (keys %{$hash_omniscient->{'level3'} } ){
				 	  							if($tag_l3 =~"utr"){
				 	  								if( exists_keys($hash_omniscient,('level3', $tag_l3, $id_l2)) ){
							 	  						foreach my $l3_feature (@{$hash_omniscient->{'level3'}{$tag_l3}{$id_l2} } ){
							 	  							if($l3_feature->_tag_value('ID') eq $UTR_location->[0][0]){
							 	  								print "UTR location modified: = ".$l3_feature->gff_string."\nnew location:".$UTRexp_location->[1]." ".$UTRexp_location->[2]."\n" if ($verbose >= 3);
							 	  								$l3_feature->start($UTRexp_location->[1]);
							 	  								$l3_feature->end($UTRexp_location->[2]);
							 	  								last;
							 	  							}
							 	  						}
							 	  					}
							 	  				}
							 	  			}
				 	  					}
				 	  					elsif($overlap and $overlap eq "perfect"){ #An UTR that match perfectly already exists !
				 	  						$create_utr=undef;
				 	  					}
				 	  				}

				 	  				if($create_utr){
				 	  					push @{$list_utr_to_create}, $new_location;
				 	  				}
			 					}
			 				}
			 				else{print "UTR check step. How is it possible ? We have an UTR in the file but none is expected according to the described exons.\nLevel2 studied:".$id_l2."\n";exit;}
		 				}
	 					else{
	 						if($list_location_UTR_expected){ 
	 							$list_utr_to_create=$list_location_UTR_expected;# no UTR exists, we have to create all of them
	 						} 
 						}
 						print "list_utr_to_create: ".Dumper($list_utr_to_create) if ($verbose >= 3); 

						# NOW CREATE UTR IF NECESSARY
						my @cds_sorted = sort {$a->[1] <=> $b->[1]} @{$list_location_CDS};

						my $extremLeftCDS = $cds_sorted[0]->[1];
						my $extremRightCDS = $cds_sorted[$#cds_sorted]->[2];

						if($list_utr_to_create){

					 	  	foreach my $location (@{$list_utr_to_create}){
					 	  		$resume_case++;
					 	  		print "_check_utrs Create one UTR !\n" if ($verbose >= 2);

					 	  		my $feature_utr = clone($feature_example);#create a copy of a random feature l3;
								$feature_utr->start($location->[1]);
					 	  		$feature_utr->end($location->[2]);
					 	  		$feature_utr->frame(".");

					 	  		#HANDLE primary tag
					 	  		my $primary_tag = "UTR";
					 	  		if($location->[2] < $extremLeftCDS){
					 	  			if($feature_utr->strand == 1){
					 	  				$primary_tag = "five_prime_UTR";
					 	  			}
					 	  			else{
					 	  				$primary_tag = "three_prime_UTR";
					 	  			}
					 	  		}
					 	  		elsif($location->[1] > $extremRightCDS){
					 	  			if($feature_utr->strand == 1){
					 	  				$primary_tag = "three_prime_UTR";
					 	  			}
					 	  			else{
					 	  				$primary_tag = "five_prime_UTR";
					 	  			}
					 	  		}

					 	  		$feature_utr->primary_tag($primary_tag);


					 	  		my $uID = _check_uniq_id($hash_omniscient, $miscCount, $uniqID, $uniqIDtoType, $feature_utr);
					 	  		create_or_replace_tag($feature_utr, 'ID', $uID); # remove parent ID because, none.
								#save new feature L2
								push (@{$hash_omniscient->{"level3"}{lc($primary_tag)}{$id_l2}}, $feature_utr);
					 	  	}
					 	}
				 	}

				 	#keep track of l2 checked (as we loop over L3, we meet several time the same l2)
		 	  		$checked{$id_l2}++;
	 	  		}
 	  		}
 	  	}
 	}
 	print "We created $resume_case UTRs that were missing\n" if($verbose >= 1 and $resume_case);
 	print "We modified $resume_case2 UTRs positions that were wrong\n" if($verbose >= 1 and $resume_case2);
}

# @Purpose: Will merge a list of "location" (tuple of integer), and another list of location. If two location overlap or are adjacent, only one location will be kept that represent the most extrem values
# @input: 3 =>  list of 3 values([[S,X,Y][S,Z,W]] or [[[S],X,Y]]),  list of integer tuple, verbose option for debug
# @output: list of list 
sub _manage_location{
	my ($locationRefList, $locationTargetList, $method, $verbose) = @_;

	_printSurrounded("Enter",25,"+","\n\n") if ($verbose >= 4); 
	print "Enter Ref: ".Dumper($locationRefList)."\nEnter Target: ".Dumper($locationTargetList) if ($verbose >= 4); 

	my @new_location_list; #new location list that will be returned once filled
	
	if ($locationTargetList){ #check number of location -> List not empty
		
		my $check_list=1;
		my $nb_ref = $#{$locationRefList};

		while($check_list){
			$check_list=undef;
			
			my %tupleChecked=();
			my %locationSaved=();
			foreach my $location_from_ref_list (@{$locationRefList}){

				if($verbose >= 4){print "\n==Location REF:".$location_from_ref_list->[1]." ".$location_from_ref_list->[2]."==\n";}
				my $location_skipped=undef; #keep track for later. If it has never been further, we will not save it later, because it will cause duplicates
				
				foreach my $location_from_target_list (@{$locationTargetList}){
					if($verbose >= 4){print "location_from_ref_list:".Dumper($location_from_ref_list)."\nAGAINST\nlocation_from_target_list:".Dumper($location_from_target_list)."\n";}

					#skip case test already done
					my @tuple = sort {$a->[1] <=> $b->[1]} ($location_from_ref_list,$location_from_target_list);
					my @tuple_ok = (@{$tuple[0]->[0]},$tuple[0]->[1],$tuple[0]->[2], @{$tuple[1]->[0]},$tuple[1]->[1],$tuple[1]->[2]);
					my $tupleString = join(' ', @tuple_ok);

					if(exists_keys (\%tupleChecked,($tupleString) ) ) {
						if($verbose >= 4){print "skip tested: $tupleString\n";}
						$location_skipped=1;
						next;
					}

					# FOLLOWING UP OF REF-TARGET / TARGET-REF test
					$tupleChecked{$tupleString}++;  # track tuple checked when we conpare a lsit aginst itself. Because a test A<=>B is equal to B<=>A
	
					#check if we modify the location or not
					my ($new_location, $overlap);
					if($method eq 'adjacent'){
					 ($new_location, $overlap) = _manage_location_lowLevel_adjacent($location_from_ref_list, $location_from_target_list);
					}else{
					 ($new_location, $overlap) = _manage_location_lowLevel_overlap($location_from_ref_list, $location_from_target_list);
					}

					my $loc = join(' ', @$new_location);
					if(! exists_keys (\%locationSaved,($loc) ) ){ #If location already saved--- we skip it
						if($verbose >= 4){print "          push1:".$new_location->[1]." ".$new_location->[2]."\n\n";}
						
						#TO PUSH THE TARGET, MODIFIED OR INTACT
						push @new_location_list, [@$new_location];
						
						$locationSaved{$loc}++;

						# FOLLOWING UP OF NewLoc-NewLoc test
						my @NewLoc_tuple = sort {$a->[1] <=> $b->[1]} ($new_location,$new_location);
						my @NewLoc_tuple_ok = (@{$NewLoc_tuple[0]->[0]},$NewLoc_tuple[0]->[1],$NewLoc_tuple[0]->[2], @{$NewLoc_tuple[1]->[0]},$NewLoc_tuple[1]->[1],$NewLoc_tuple[1]->[2]);
						my $NewLoctString = join(' ', @NewLoc_tuple_ok);
						$tupleChecked{$NewLoctString}++;
					}

					# FOLLOWING UP OF Target-Target test
					my @TargetTarget_tuple = sort {$a->[1] <=> $b->[1]} ($location_from_target_list,$location_from_target_list);
					my @TargetTarget_tuple_ok = (@{$TargetTarget_tuple[0]->[0]},$TargetTarget_tuple[0]->[1],$TargetTarget_tuple[0]->[2], @{$TargetTarget_tuple[1]->[0]},$TargetTarget_tuple[1]->[1],$TargetTarget_tuple[1]->[2]);
					my $TargetTargetString = join(' ', @TargetTarget_tuple_ok);
					$tupleChecked{$TargetTargetString}++;

					if( ($new_location->[1] !=  $location_from_target_list->[1] or $new_location->[2] !=  $location_from_target_list->[2]) ){ # location has been modified or not modifier but overlap (It means overlap completly ... it is included in)
						$check_list=1;
						if($verbose >= 4){print "LOCATION MODIFIED: ".$location_from_target_list->[2]." ".$new_location->[2]."\n";}
					}
					elsif($overlap){#position not modified but overlap (A is completely included in B); Need to keep track of it to not save the position when we are out of the loop
						$location_skipped=1; 
					}
				}

				#TO PUSH THE REF
				if(! $check_list and ! $location_skipped){ #No modification done
					my $loc = join(' ', @$location_from_ref_list);

					if(! exists_keys (\%locationSaved,($loc) ) ){
						if($verbose >= 4){print " push LocationREF = add new value !!\n";}
						push @new_location_list, [@{$location_from_ref_list}];
					}
				}
			}

			#modification done we have to re-check all location between each other
			if($check_list){ 

				if (scalar @new_location_list == 1){
					$check_list = undef; #Do not check the list if it is size one ! Because we will be stuck => case test againt himself is skipped !
				}
				elsif($nb_ref > 0){ #Case where the target list as input was more than 1 element, we have to check all against all  due to case like Ref=> -------    ---------
					$locationTargetList = [@new_location_list];																						# Tearget ------------
					$locationRefList = [@new_location_list];
					$nb_ref = $#{$locationRefList}; #reinitialise the value
					if($verbose >= 4){print "Location in memory:".@new_location_list." ".Dumper(\@new_location_list)." NNNNNNNNNNNNNNNNNNNNNNNow check aginst itself !\n";}
					@new_location_list=();
					%tupleChecked=();
				}
				else{$check_list=undef;}					
			}
		}
	}
	else{#check number of location -> none
		_printSurrounded("Return",25,"-","\n\n") if ($verbose >= 4);
		if($verbose >= 4){print "returnA: ".Dumper($locationRefList)."\n\n\n";}
		return \@{$locationRefList};
	}
	_printSurrounded("Return",25,"-","\n\n") if ($verbose >= 4); 
	if($verbose >= 4){print "returnB: ".Dumper(\@new_location_list)."\n\n\n";}
	return \@new_location_list;
}

#	===================== location1
#		===================== location2    
#   ========================= <= New location2 returned
# @Purpose: Modify the location2 if it overlap the location1 by keeping the extrem values. Return the location2 intact if no overlap. /!\ The locations are merged if they are contigu
# @input: 2 =>  integer tuple [[ID],X,Y],  list of integer tuple
# @output: 2 => ref of a list of 2 element, boolean
sub _manage_location_lowLevel_adjacent{
	my ($location, $location2) = @_;

	my $new_location = [@{$location2}];
	my $overlap=undef;

	if ( ($location2->[1] <= $location->[2]+1) and ($location2->[2]+1 >= $location->[1]) ){ #it overlaps or are consecutive
		
		#Manage Id to avoid same IDs
		my %params = map { $_ => 1 } @{$new_location->[0]};
		foreach my $id ( @{$location->[0]}){
			if(! exists($params{$id})){ 
				push @{$new_location->[0]}, $id ; #append from the end the list of ID
			}
		}
		$overlap=1;

		if($location2->[1] > $location->[1]){
			$new_location->[1]=$location->[1];
		}
		if($location->[2] > $location2->[2]){
			$new_location->[2]=$location->[2];
		}
	}
	return $new_location, $overlap;
}

#	===================== location1
#		===================== location2    
#   ========================= <= New location2 returned
# @Purpose: Modify the location2 if it overlap the location1 by keeping the extrem values. Return the location2 intact if no overlap. /!\ We append the ID list by the end (as push) when there is an overlap
# @input: 2 =>  integer tuple [[ID],X,Y],  list of integer tuple
# @output: 2 => ref of a list of 2 element, boolean
sub _manage_location_lowLevel_overlap{
	my ($location, $location2) = @_;

	my $new_location = [@{$location2}];
	my $overlap=undef;

	if ( ($location2->[1] <= $location->[2]) and ($location2->[2] >= $location->[1]) ){ #it overlaps or are consecutive

		#Manage Id to avoid same IDs
		my %params = map { $_ => 1 } @{$new_location->[0]};
		foreach my $id ( @{$location->[0]}){
			if(! exists($params{$id})){ 
				push  @{$new_location->[0]}, $id ; #append from the end the list of ID
			}
		}

		$overlap=1;

		if($location2->[1] > $location->[1]){
			$new_location->[1]=$location->[1];
		}
		if($location->[2] > $location2->[2]){
			$new_location->[2]=$location->[2];
		}
	}

	return $new_location, $overlap;
}


#	================= 		  location1 (cds)
#		===================== location2 (exon)   
#                    ======== <= New location2 returned
sub _manage_location_lowLevel_inversed{
	my ($location, $location2, $verbose) = @_;
	
	print "_manage_location_lowLevel_inversed\n" if($verbose >= 3);

	my $new_location = [@{$location2}];
	my $overlap=undef;

	if ( ($location2->[1] == $location->[1]) and ($location2->[2] == $location->[2]) ){ #it overlaps perfectly
		return $new_location, "perfect";
	}

	if ( ($location2->[1] <= $location->[2]) and ($location2->[2] >= $location->[1]) ){ #it overlaps

		$overlap=1;

		if($location2->[1] < $location->[1]){
			$new_location->[2] = $location->[1]-1;
			print "case1\n" if($verbose >= 3);
		}
		if($location->[2] < $location2->[2]){
			$new_location->[1] = $location->[2]+1;
			print "case2\n" if($verbose >= 3);
		}
	}
	return $new_location, $overlap;
}

#============================================================================================================
#Explanation: Case where part of the locus BBBBBB has been seen before to meet the its Parent feature (see below) = a parent feature ID has been created on the fly during the parsing.
#			  We now need to remove the wrong Parent ID and link them to the correct one.
#seq1	maker	CDS	561401	561519	.	+	2	ID=CLUHART00000006146:cds;locus_tag=AAAAA
#seq1	maker	UTR	337818	337914	.	+	.	ID=CLUHART00000008717:five_prime_utr;locus_tag=BBBBBB
#seq1	maker	UTR	343034	343277	.	+	.	ID=CLUHART00000008717:three_prime_utr;locus_tag=BBBBBB
#seq1	maker	CDS	564171	564235	.	+	0	ID=CLUHART00000006146:cds;locus_tag=AAAAA
#...
#seq1	maker	gene	337818	343277	.	+	.	ID=CLUHARG00000005458;locus_tag=BBBBBB
#
# HIS<=>Hash InfoSequential
#
# => Improve something ?????? 
sub _cleanSequentialIncase{
	my ($infoSequential, $locusTAGuniq, $verbose) = @_;
	my $resume_case=undef;

	foreach my $locusNameHIS (keys %{$infoSequential} ){

	 	if(exists_keys($locusTAGuniq,('level1', $locusNameHIS))){
	 		
	 		my $locusNameUniq = $locusTAGuniq->{'level1'}{$locusNameHIS};
	 		if($locusNameHIS ne $locusNameUniq ){
	 			$resume_case++;
	 			
	 			# The locusNameUniq already exists, we have to fill it with the part of inforamtion missing that is contained in$infoSequential->{$locusNameHIS}
	 			if(exists_keys ($infoSequential,($locusNameUniq) ) ){
	 				
	 				foreach my $bucket (keys %{$infoSequential->{$locusNameHIS}} ){
	 					if ($bucket eq 'level1'){next;}
	 					
	 					my $prefix= lc(PREFIXL2); #when a l2 start with this prefix it means we created the l2 on the fly (the real l2 if exists, had not been met yet)
	 					if($bucket =~ /^$prefix/i){
	 						my $idok=undef;
	 						foreach my $feature ( @{$locusTAGuniq->{'level2'}{ $locusNameUniq}}){	 							
	 							if(lc($feature->_tag_value('ID')) !~ /^$prefix/i){
	 								$idok = lc( $feature->_tag_value('ID') ); # @{$locusTAGuniq->{'level2'}{ $locusNameUniq }}[$cpt] is the first l2 feature that has been realy met 
	 								last;									  # We make the assumption that the pieces of the locus that were lost before to describe its real l2 is part of the first real l2 met.
	 																	      # ====================================================================================================================================						
	 							}
	 						}

	 				 		if(exists_keys ($infoSequential,($locusNameUniq, $idok) ) ){
	 				 			
	 				 			foreach my $level (keys %{$infoSequential->{$locusNameHIS}{$bucket}} ){
	 				 				push @{$infoSequential->{$locusNameUniq}{$idok}{$level}}, @{$infoSequential->{$locusNameHIS}{$bucket}{$level}};
	 				 			}
	 				 			delete $infoSequential->{$locusNameHIS}{$bucket};
	 				 			if(! %{$infoSequential->{$locusNameHIS}}){delete $infoSequential->{$locusNameHIS};} # remove because nothing linked to it anymore
	 				 		}
	 				 		else{
	 				 		  $infoSequential->{$locusNameUniq}{$idok} = delete $infoSequential->{$locusNameHIS}{$bucket}; #delete the lod key but transfer the data to a new key
	 				 		}
	 				 	}
	 			 	}
	 			}
	 			else{ # The locusNameUniq didn't exists, we can directly shift the old locusNameHIS by the locusNameUniq
	 				$infoSequential->{$locusNameUniq} = delete $infoSequential->{$locusNameHIS}; # link to the first l2 #delete the lod key but transfer the data to a new key
	 			}
			}
	 	}
	}
	print "We found $resume_case cases where part of the locus data where defined earlier in the file.\n" if($verbose >= 1 and $resume_case);
}

#
#
# All Level1 feature are for sure in omniscient, we have only the ID in sequential, other level feature are in sequential only if no parent has been found
sub _check_sequential{ # Goes through from L3 to l1
 	my ($infoSequential, $omniscient, $miscCount, $uniqID, $uniqIDtoType, $locusTAGuniq, $mRNAGeneLink, $verbose) = @_;
 	my $resume_case=undef;

 	_cleanSequentialIncase($infoSequential, $locusTAGuniq, $verbose); # PART OF LOCUS LOST BEFORE TO MEET IT L2 or L1 ... we catch them and re-link everythong as it should be

 	foreach my $locusNameHIS (keys %{$infoSequential} ){ #comon tag was l1 id wheb no real comon tag present

 		foreach my $bucket (keys %{$infoSequential->{$locusNameHIS} } ){ #bucket = level1 or Id L2
 			print "\nlocusNameHIS $locusNameHIS bucket $bucket\n\n" if ($verbose >= 3);
 			
 			if ($bucket eq 'level1'){next;} #skip case level1 - structure of the hash different

 			my $must_create_l2=undef;
 			my $feature_l2 = undef;

			#Bucket is an uniq ID created during the reading process. So it can be used as uniq ID.
 			if(! exists_keys($infoSequential,($locusNameHIS, $bucket, 'level3') ) ){
 				
 				# Link the l2 to the L1 feature
				$feature_l2=$infoSequential->{$locusNameHIS}{$bucket}{'level2'};
				print "level2 in sequential doenst have L3 feature associated in sequential - $locusNameHIS $bucket! ".$feature_l2->gff_string."\n" if ($verbose >= 3);

				if(! exists($mRNAGeneLink->{lc($bucket)}) ){
					# We add it to omniscient and to mRNAGeneLink
					print "level2 does not exits in mRNAGeneLink(omniscient) !".$feature_l2->gff_string."\n" if ($verbose >= 3);
					push (@{$omniscient->{"level2"}{lc($feature_l2->primary_tag)}{lc($feature_l2->_tag_value('Parent'))} }, $feature_l2);
					$mRNAGeneLink->{lc($feature_l2->_tag_value('ID'))} = $feature_l2->_tag_value('Parent');
				}
				if( ! exists_keys($omniscient,('level3', "exon", lc($feature_l2->_tag_value("ID")))) ){ #check if an exon exist in the omniscient
					if ($createL3forL2orphan){ # create the exon missing if option agreed
		 				print "create single level3 exon feature  !\n" if($verbose >= 2);
		 				my $feature_l3 = clone($feature_l2);#create a copy of the l2 feature;
		 				$feature_l3->primary_tag('exon');
		 				create_or_replace_tag($feature_l3,'Parent', $feature_l3->_tag_value('ID')); # change parentID
						$feature_l3->remove_tag('ID');
						# create ID
						my $id =  _create_ID($miscCount, $uniqID, $uniqIDtoType, 'exon', "exon-1", 'nbis_NEW');
						create_or_replace_tag($feature_l3,'ID', $id); # change ID
		    			#my $id = _check_uniq_id($omniscient, $miscCount, $uniqID, $uniqIDtoType, $feature_l3);
		    			push (@{$omniscient->{"level3"}{lc($feature_l3->primary_tag)}{lc($feature_l3->_tag_value('Parent'))} }, $feature_l3);
					}
				}
 				#warn "Not normal, we have feature L2  without L3 feature associated.\n"; #We cannot guess the structure except if it is prokaryote or single exon in eucaryote... should we improve that ?
 			}
			else{
 				foreach my $feature_L3 (@{$infoSequential->{$locusNameHIS}{$bucket}{'level3'}} ){

 					if(! exists_keys($infoSequential,($locusNameHIS, $bucket,'level2'))  ){
 						print "_check_sequential level2 does not exits in sequential !\n" if($verbose >= 2);

 						#take L2 from omniscient if already exits
 						if(exists($mRNAGeneLink->{lc($bucket)}) ){

 							my $l1_id = $mRNAGeneLink->{lc($bucket)};
 							foreach my $tag_l2 (keys %{$omniscient->{'level2'}}){
 								if(exists_keys($omniscient, ('level2', $tag_l2, lc($l1_id) ) ) ){
		 							foreach my $featureL2 (@{$omniscient->{'level2'}{$tag_l2}{lc($l1_id)}}){
		 								if(lc($featureL2->_tag_value('ID')) eq $bucket){
		 									print "_check_sequential level2 exits in omniscient !\n" if($verbose >= 2);
		 									$feature_l2 = $featureL2;
		 									last;
		 								}
		 							}
		 							if($feature_l2){last;}
		 						}
	 						}
 							$infoSequential->{$locusNameHIS}{$bucket}{'level2'} = $feature_l2;
						}
 						else{#create l2
 							print "create level2  !\n" if($verbose >= 2);
	 						$must_create_l2=1;
	 						$feature_l2 = clone($infoSequential->{$locusNameHIS}{$bucket}{'level3'}[0]);#create a copy of the first mRNA feature;
							
							#manage primary tag
							my $primary_tag_l2='RNA';
							foreach my $feature_L3 (@{$infoSequential->{$locusNameHIS}{$bucket}{'level3'}} ){

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
							 	if( exists_keys($infoSequential,($locusNameHIS,'level1'))  ){ # parent ID exists in infoSequential
	 								$parentID = lc($infoSequential->{$locusNameHIS}{'level1'}); # PArentID it correct case ???
	 								print "_check_sequential Parent IDtaken from infoSequential\n" if ($verbose >= 3);
	 							}
								else{
									my $IDgoodCast = _id_exists_in_l1_omniscient($omniscient, $locusNameHIS);
									if($IDgoodCast){
											$parentID = $IDgoodCast;
											print "_check_sequential Parent IDtaken from omniscient\n" if ($verbose >= 3);
									}

						 			if( ! $parentID ){ #In that case level1 feature doesn't exists in $infoSequential and in $omniscient. I will be created by the method check_gene_link_to_mrna 
						 				#my	($miscCount, $uniqID, $primary_tag, $id, $prefix)=@_;
						 				$parentID =  _create_ID($miscCount, $uniqID, $uniqIDtoType, 'gene', "gene-1", 'nbis_NEW');
						 				print "_check_sequential Parent IDtaken created\n" if ($verbose >= 3);
						 				$infoSequential->{$locusNameHIS}{'level1'}=$parentID;
						 			}
						 		}
						 		print "_check_sequential Parent ID created for level2 = $parentID\n" if ($verbose >= 2);
					 			create_or_replace_tag($feature_l2,'Parent', $parentID ); # change parentID 

					 		print "push-omniscient: level2 || ".lc($primary_tag_l2)." || ".lc($parentID)." == ".$feature_l2->gff_string."\n" if ($verbose >= 2);
					 		push (@{$omniscient->{"level2"}{lc($primary_tag_l2)}{lc($parentID)}}, $feature_l2);
					 		$mRNAGeneLink->{lc($bucket)} = $parentID; # Always need to keep track about l2->l1, else the method check_l3_link_to_l2 will recreate a l1 thinking this relationship is not fill
					 		$infoSequential->{$locusNameHIS}{$bucket}{'level2'} = $feature_l2;
					 	}
 					}
					else{
						
						#MUST push L2 in omniscient if absent !
						$feature_l2=$infoSequential->{$locusNameHIS}{$bucket}{'level2'};
						print "level2 exits in sequential - $locusNameHIS $bucket! ".$feature_l2->gff_string."\n" if ($verbose >= 3);

						if(! exists($mRNAGeneLink->{$bucket}) ){
							print "level2 does not exits in mRNAGeneLink(omniscient) !".$feature_l2->gff_string."\n" if ($verbose >= 3);
							push (@{$omniscient->{"level2"}{lc($feature_l2->primary_tag)}{lc($feature_l2->_tag_value('Parent'))} }, $feature_l2);
							$mRNAGeneLink->{lc($feature_l2->_tag_value('ID'))} = $feature_l2->_tag_value('Parent');
						}
					}					
 					
 					my $primary_tag_L3 =  lc($feature_L3->primary_tag);
 					create_or_replace_tag($feature_L3,'Parent', $feature_l2->_tag_value('ID')); #modify ID to replace by parent value

 					print "push-omniscient: level3 || ".$primary_tag_L3." || ".$bucket." == ".$feature_L3->gff_string."\n" if ($verbose >= 2);
 					push (@{$omniscient->{"level3"}{$primary_tag_L3}{$bucket}}, $feature_L3);
 				}
 				$resume_case++;
 			}

 			if($must_create_l2){
 				check_level2_positions($omniscient, $feature_l2);
 			}
 		}
 		#LEVEL 1 IS taking care later
 	}
 	print "We found $resume_case sequential cases\n" if($verbose >= 1 and $resume_case);
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


# Check the start and end of level1 feature based on all features level2;
sub _check_all_level1_positions {
	my ($hash_omniscient, $verbose)=@_;
	my $resume_case=undef;

	foreach my $tag_l1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
		foreach my $id_l1 ( keys %{$hash_omniscient->{'level1'}{$tag_l1}} ) { #sort by position
			
			my $level1_feature = $hash_omniscient->{'level1'}{$tag_l1}{$id_l1};

			$resume_case++ if(check_level1_positions($hash_omniscient, $level1_feature, $verbose));
		}
	}
	print "We fixed $resume_case wrong level1 position cases\n" if($verbose >= 1 and $resume_case);
}

# Purpose: review all the feature L2 to adjust their start and stop according to the extrem start and stop from L3 sub features.
sub _check_all_level2_positions{
	my ($hash_omniscient, $verbose)=@_;
	my $resume_case=undef;

	foreach my $tag_l1 (keys %{$hash_omniscient->{'level1'}}){ # primary_tag_key_level1 = gene or repeat etc...
		foreach my $id_l1 ( keys %{$hash_omniscient->{'level1'}{$tag_l1}} ) { #sort by position
			
			foreach my $tag_level2 (keys %{$hash_omniscient->{'level2'}}){ # primary_tag_key_level2 = mrna or mirna or ncrna or trna etc...
    			if ( exists_keys ($hash_omniscient, ('level2', $tag_level2, $id_l1) ) ){
					
					foreach my $mRNA_feature ( @{$hash_omniscient->{'level2'}{$tag_level2}{$id_l1}}){ 
						my $level2_ID = lc($mRNA_feature->_tag_value('ID'));
						my @feature_list=();	
						foreach my $primary_tag_l3 ( keys %{$hash_omniscient->{'level3'}}){ # primary_tag_l3 = cds or exon or start_codon or utr etc...
							
							if ( exists_keys( $hash_omniscient, ('level3', $primary_tag_l3, $level2_ID) ) ){
								push @feature_list, @{$hash_omniscient->{'level3'}{$primary_tag_l3}{$level2_ID}};
							}
						}
						if(scalar(@feature_list) > 0){ #could be emtpy like in match match_part features
							 $resume_case++ if(check_mrna_positions($mRNA_feature, \@feature_list, $verbose));
						}
					}
    			}
    		}
		}
	}

	print "We fixed $resume_case wrong level2 position cases\n" if($verbose >= 1 and $resume_case);
}

# Check the start and end of mRNA based a list of feature like list of exon;
sub check_mrna_positions{
  my ($mRNA_feature, $exon_list, $verbose)=@_;
  if(! $verbose){$verbose=0;}
  my $result=undef;

  my @exon_list_sorted = sort {$a->start <=> $b->start} @{$exon_list};
  my $exonStart=$exon_list_sorted[0]->start;

  @exon_list_sorted = sort {$a->end <=> $b->end} @exon_list_sorted;
  my $exonEnd=$exon_list_sorted[$#exon_list_sorted]->end;

  #check start
  if ($mRNA_feature->start != $exonStart){
  	print "We modified the L2 LEFT extremity for the sanity the biological data!\n" if($verbose >= 3 );
    $mRNA_feature->start($exonStart);
    $result=1;
  }
  #check stop
  if($mRNA_feature->end != $exonEnd){
  	print "We modified the L2 RIGHT extremity for the sanity the biological data!\n" if($verbose >= 3);
    $mRNA_feature->end($exonEnd);
    $result=1;
  }

  return $result;
}

# L1: LocusID->level->typeFeature->ID->[ID,start,end]
# LocusID->level->typeFeature->Parent->[ID,start,end]
# @Purpose: When two feature overlap at level3, and are the same type level 2 they have to be merged under the same level 1 feature.
# @input: 2 =>  hash,  integer for verbosity
# @output: 0 
sub _check_overlap_name_diff{
	my ($omniscient, $mRNAGeneLink, $verbose) = @_;
	my $resume_case=undef;

	my $sortBySeq = _sort_by_seq($omniscient);
	my %alreadyChecked;

	foreach my $locusID ( keys %{$sortBySeq}){ # tag_l1 = gene or repeat etc...
		if ( exists_keys( $sortBySeq, ( $locusID, 'level1') ) ){
			foreach my $tag_l1 ( keys %{$sortBySeq->{$locusID}{'level1'}} ) { 

				my $to_check = clone($sortBySeq->{$locusID}{'level1'}{$tag_l1});

				# Go through location from left to right ### !!
				foreach my $id_l1 ( sort {$sortBySeq->{$locusID}{'level1'}{$tag_l1}{$a}[1] <=> $sortBySeq->{$locusID}{'level1'}{$tag_l1}{$b}[1] } keys %{$sortBySeq->{$locusID}{'level1'}{$tag_l1}} ) {

					if(! exists($alreadyChecked{$id_l1})){

						#remove itself 
						delete $to_check->{$id_l1};
						my $location = $sortBySeq->{$locusID}{'level1'}{$tag_l1}{$id_l1}; # This location will be updated on the fly

						# Go through location from left to right ### !!
						foreach my $id2_l1 ( sort {$to_check->{$a}[1] <=> $to_check->{$b}[1] } keys %{$to_check} ) {

							my $location_to_check = $to_check->{$id2_l1};

							#If location_to_check start if over the end of the reference location, we stop
							if($location_to_check->[1] > $sortBySeq->{$locusID}{'level1'}{$tag_l1}{$id_l1}[2]) {last;} 

							my ($location, $overlap) = location_overlap_update($location, $location_to_check); # location is updated on the fly, and the newly modified location is the one that will be used at the next loop
							
							# Let's check at Gene LEVEL
							if($overlap){

								#let's check at CDS level
								if(check_gene_overlap_at_CDSthenEXON($omniscient, $omniscient , $id_l1, $id2_l1)){ #If contains CDS it has to overlap at CDS level to be merged, otherwise any type of feature level3 overlaping is sufficient to decide to merge the level1 together
									#they overlap in the CDS we should give them the same name
									$resume_case++;
									$alreadyChecked{$id_l1}++;

									print "$id_l1 and  $id2_l1 same locus. We merge them together. Below the corresponding feature groups in their whole.\n" if ($verbose >= 3);
									print_omniscient_from_level1_id_list($omniscient, [$id_l1,$id2_l1], $fh_error ) if ($verbose >= 3);
									# remove the level1 of the ovelaping one
									delete $omniscient->{'level1'}{$tag_l1}{$id2_l1};
									# remove the level2 to level1 link stored into the mRNAGeneLink hash. The new links will be added just later after the check to see if we keep the level2 feature or not (we remove it when identical)
									foreach my $l2_type (%{$omniscient->{'level2'}}){
										if(exists_keys($omniscient,('level2', $l2_type, $id2_l1))){
											foreach my $feature_l2 (@{$omniscient->{'level2'}{$l2_type}{$id2_l1}}){
												delete $mRNAGeneLink->{lc($feature_l2->_tag_value('ID'))};
											}
										}
									}

									# Let's change the parent of all the L2 features 
									foreach my $l2_type (%{$omniscient->{'level2'}} ){

										if(exists_keys($omniscient,('level2', $l2_type, $id2_l1))){
											###############################
											# REMOVE THE IDENTICAL ISOFORMS

											# first list uniqs
											my $list_of_uniqs  = keep_only_uniq_from_list2($omniscient, $omniscient->{'level2'}{$l2_type}{$id_l1}, $omniscient->{'level2'}{$l2_type}{$id2_l1}, $verbose); # remove if identical l2 exists
													

											#Now manage the rest
											foreach my $feature_l2 (@{$list_of_uniqs}){

												create_or_replace_tag($feature_l2,'Parent', $sortBySeq->{$locusID}{'level1'}{$tag_l1}{$id_l1}[0]); #change the parent
												# Add the corrected feature to its new L2 bucket
												push (@{$omniscient->{'level2'}{$l2_type}{$id_l1}}, $feature_l2);
												
												# Attach the new parent into the mRNAGeneLink hash
												$mRNAGeneLink->{lc($feature_l2->_tag_value('ID'))}=$feature_l2->_tag_value('Parent');

											}
											# remove the old l2 key
											delete $omniscient->{'level2'}{$l2_type}{$id2_l1};
										}
									}
									check_level1_positions($omniscient, $omniscient->{'level1'}{$tag_l1}{$id_l1}, 0);

									if($omniscient->{'level1'}{$tag_l1}{$id_l1}->end > $sortBySeq->{$locusID}{'level1'}{$tag_l1}{$id_l1}[2]){
										$sortBySeq->{$locusID}{'level1'}{$tag_l1}{$id_l1}[2] = $omniscient->{'level1'}{$tag_l1}{$id_l1}->end;
										print "This one Now !!\n";exit;
									}
								}
							}
						}
					}
				}
	 		}
	 	}
	}
	print "We fixed $resume_case case where feature has been merged within the same locus\n" if($verbose >= 1 and $resume_case);
}


# @Purpose: When too feature l2 isoform are identical, we remove one
# @input: 2 =>  hash,  integer for verbosity
# @output: 0 
sub _check_identical_isoforms{
	my ($omniscient, $mRNAGeneLink, $verbose) = @_;
	my $resume_case=undef;

	# Go through oall l2 feature
	foreach my $l2_type (keys %{$omniscient->{'level2'}}){
		foreach my $id2_l1 (keys %{$omniscient->{'level2'}{$l2_type}}){
			# If more than 1 related to level1
			
			if(exists_keys($omniscient,('level2', $l2_type, $id2_l1)) and scalar @{$omniscient->{'level2'}{$l2_type}{$id2_l1}} > 1){ # more than one l2 feature of that type

				my @L2_list_to_remove;
				my %checked;
				foreach my $feature2 (sort {$b cmp $a} @{$omniscient->{'level2'}{$l2_type}{$id2_l1}}){			
					$checked{lc($feature2->_tag_value('ID'))}{lc($feature2->_tag_value('ID'))}++;

					my $keep = 1;
					foreach my $feature1 (sort {$b cmp $a} @{$omniscient->{'level2'}{$l2_type}{$id2_l1}}){

						# If not itself and not already checked (A -> B is the same as B -> A), and A or B already removed and must now be skiped (skipme key)
						if( (! exists_keys(\%checked, (lc($feature2->_tag_value('ID')), "skipme"))) and (! exists_keys(\%checked, (lc($feature1->_tag_value('ID')), "skipme"))) and ! exists_keys(\%checked, ( lc($feature2->_tag_value('ID')), lc($feature1->_tag_value('ID')) ) ) ){ #
							$checked{lc($feature2->_tag_value('ID'))}{lc($feature1->_tag_value('ID'))}++;
							$checked{lc($feature1->_tag_value('ID'))}{lc($feature2->_tag_value('ID'))}++;

							#check their position are identical
							if($feature1->start().$feature1->end() eq $feature2->start().$feature2->end()){

								#Check their subfeature are  identicals
								if(l2_identical($omniscient, $feature1, $feature2, $verbose )){
									$keep = undef; 
									last;
								}
							}
						}
					}
					# We dont keep the l2 feature so we have to remove all related features and itself
					if(! $keep){ 
						print "Lets remove isoform ".$feature2->_tag_value('ID')."\n" if ($verbose == 2);
						$checked{lc($feature2->_tag_value('ID'))}{"skipme"}++;# will be removed later do not check anymore this one
						
						foreach my $tag (keys %{$omniscient->{'level3'}}){
							if(exists_keys($omniscient, ('level3', $tag, lc($feature2->_tag_value('ID'))))){
								delete $omniscient->{'level3'}{$tag}{lc($feature2->_tag_value('ID'))};
							}
						}
						#Hast to be removed once we finished to go through the l2 list
						push(@L2_list_to_remove,lc($feature2->_tag_value('ID')));
						delete $mRNAGeneLink->{lc($feature2->_tag_value('ID'))};
					}
				}

				#L2 has to be removed from List
				my @newL2List;	
				foreach my $feature ( @{$omniscient->{'level2'}{$l2_type}{$id2_l1}} ){
					my $keep = 1;
					foreach my $id_l2 (@L2_list_to_remove){
						if( lc($feature->_tag_value('ID')) eq lc($id_l2) ){
							$keep = undef;
						}
					}
					if($keep){
						push (@newL2List,$feature)
					}
				}
				@{$omniscient->{'level2'}{$l2_type}{$id2_l1}}=@newL2List;

			}
		}
	}
}


#
# Sort by locusID !!!!
# L1 => LocusID->level->typeFeature->ID =[ID,start,end]
# L2 and L3 => LocusID->level->typeFeature->Parent->ID = [ID,start,end]
#
#
sub _sort_by_seq{
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


# Actually the duplicates have been collected during the parsing process here we just print them.
sub _check_duplicates{
	my ($duplicate, $omniscient, $verbose) = @_  ;

	my $keyExist = keys %{$duplicate};
    if($keyExist){#print result
    	_printSurrounded("Achthung /\\ Attention /\\ Be carefull => Duplicates removed !\n(Same chr/contig/scaffold, same position, same ID)",75, "#");
      	
      	my $gffout= Bio::Tools::GFF->new( -fh => \*STDOUT );
      	my $info = _print_duplicates($duplicate, $omniscient, $gffout, $verbose);
    	print "$info\n" if($verbose >= 1);
    }
}

# print duplicate hash
sub _print_duplicates {
	my ($duplicate_omniscient, $hash_omniscient, $gffout, $verbose) = @_  ;

	my $string="";
	foreach my $level (keys %{$duplicate_omniscient}){ # primary_tag_key_level1 = gene or repeat etc...
		foreach my $primary_tag (keys %{$duplicate_omniscient->{$level}}){
			my $nb_by_pt=0;
			my $nb_feat_pt=0;
			foreach my $id (keys %{$duplicate_omniscient->{$level}{$primary_tag}}){
				$nb_feat_pt++;
				foreach my $feature (@{$duplicate_omniscient->{$level}{$primary_tag}{$id}}){
					$nb_by_pt++;
					$gffout->write_feature($feature) if($verbose >= 2); # print feature
				}
			}
			$string .= "There were $nb_feat_pt duplicated $primary_tag feature for a total of $nb_by_pt duplicates.";
		}
	}

	return $string;
}


# allows to add a frame to a string to print
sub _printSurrounded{
  my ($term,$size,$char,$extra) = @_;

   my $frame=$char x ($size+4);
  $frame.="\n";

  my $result = $frame;

  my @lines = split(/\n/,$term);

  	foreach my $line (@lines){
  		$result .="$char ";

  		my $sizeTerm=length($line);
	  	if ($sizeTerm > $size ){	    
		    $result .= substr($line, 0,($size));#
	 	 }
	 	else{
		    my $nbBlancBefore=int(($size-$sizeTerm) / 2);
		    my $nbBlancAfter = ($size-$sizeTerm) - $nbBlancBefore;
		    $result .= " " x $nbBlancBefore;
		    $result .= $line;
		    $result .= " " x $nbBlancAfter;	       		   
	  	}
	  	$result .= " $char\n";
	}
	$result .= "$frame";
	if($extra){$result .= "$extra";}
	print $result;
}

#GFF format guess
# Input: filename
# Output: Integer (1,2 or 3)
sub select_gff_format{
    my ($file) = @_;

    #HANDLE format
    my $format=3;
    my $problem3=undef;
    my $nbLineChecked=100; #number line to use to check the formnat
    my $cpt=0;
    
    open(my $fh, '<', $file) or die "cannot open file $file";
    {
        while(<$fh>){

        	if($_ =~ /^#/){next;} #if it is a commented line starting by # we skip it.
        	
            $cpt++;
            if($cpt > $nbLineChecked){
                    _printSurrounded("Doesn't look like a GFF file\nLet's see what the Bioperl parser can do with that...",100,"!");  
                    last;
            }
            if($_ =~ /^.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*\t(.*)/){
                if(length($1) < 1){next;}
                
                my $string = $1;
                if($string =~ /=/  and $string =~ /;/ ){ last; };

                if($string !~ /=/  and $string !~ /;/ ){  
                        $format=1;
                        #_printSurrounded("Problem detected wihtin the 9th colum of the gff file\nYou cannot have space between attributes and between tag and values.\nAll your attributes will be gathered within the GROUP tag",100,"!");  
                        last;
                }
                elsif($string !~ /=/  and $string =~ /;/ ){       
                                $format=2;
                                last;
                }
                my $c = () = $string =~ /=/g;
                my $d = () = $string =~ /\ /g;
                if($c > 1 and $d > 1  and $string !~ /;/ ){
                       $problem3=1;
                }
     	   }
        }
    }
    close($fh);
    if($problem3){
        _printSurrounded("Thre is a problem with your GFF format.\nThis format is wrong: tag=value tag=value.\nYou should have: tag=value;tag=value or tag value ; tag value\nThe best parser we can use will keep only the first attribute.",100,"!");  
    	$format=1;
    }
    return $format;
}

# We modify the attributes: group=gene_id "e_gw1.5.2.1" protein_id 335805 exonNumber 1
# in order to get : gene_id=e_gw1.5.2.1;protein_id=335805;exonNumber=1
sub _gff1_corrector{
	my ($feat)=@_;

	if($feat->has_tag('group')){
		my @attribs = $feat->get_tag_values('group');
		my $attribs = join ' ', @attribs;
		my @parsed;
		my $flag = 0; # this could be changed to a bit and just be twiddled

	    # run through each character one at a time and check it
	    # NOTE: changed to foreach loop which is more efficient in perl
	    # --jasons
	    my $previousChar=undef;
	    my $string="";
	    for my $a ( split //, $attribs ) { 
	    	$string.=$a;
	        # flag up on entering quoted text, down on leaving it
	        if( $a eq '"') { $flag = ( $flag == 0 ) ? 1:0 } #active deactive the flag
	        if ($previousChar and $previousChar eq '"' and $flag == 0){ # case we have to strip the " characters
	        	chop $string;
	        	chop $string;
	        	$string = reverse($string);
    			chop($string);
    			$string= reverse($string);
	        	push @parsed, $string;
	        	$string="";
	        }
	        elsif( $a eq " " and $flag == 0){
	        	chop $string;
	        	push @parsed, $string;
	        	$string="";
	        }
	        $previousChar = $a;
	    }

	    if($string != ""){ if($previousChar eq " "){chop $string;} push @parsed, $string;}
	    while (@parsed){
	    	my $value = pop @parsed;
	    	my $tag = pop @parsed;
	    	$feat->add_tag_value($tag, $value);
	    }
	    #remove it   
		$feat->remove_tag('group');
    } 
}

# @Purpose: Create a hash containing all the name and identifier of an ontology.
# @input: 1 =>  Object Bio::Ontology
# @output: 1 => hash containing all the name and identifier
sub create_term_and_id_hash{
    my ($self) = @_;

    my %hash_term_id;

    foreach my $term ($self->get_all_terms) {
       $hash_term_id{lc($term->name)} = lc($term->identifier);
       $hash_term_id{lc($term->identifier)} = lc($term->name);
       #print $term->name." <=> ".$term->identifier."\n";
    }
    return \%hash_term_id;
}

#Look for gff3 specific header 
#@INPUT: 1 => string (a file)
#@OUPUT: 1 => hash of the different header and their values
sub _check_header{
    my ($file) = @_;

    #HANDLE format
    my %headerInfo;
    
    #check it is a file
    if(-f $file){
	    open(my $fh, '<', $file) or die "cannot open file $file";
	    {
	        while(<$fh>){
	            if($_ !~ /^##[^#]/) {
	                 last;
	            }
	            else{
	            	my @data = split /\s/, $_ ;
	            	my $type = shift @data;

	            	if($type eq /^##gff-version/){
	            		$headerInfo{$type}=$data[0]; #1 element
	            	}
					if($type eq "##sequence-region"){
						$headerInfo{$type}=@data; # 3 elements
					}
					if($type eq "##feature-ontology"){
						$headerInfo{$type}=$data[0] #1 element
					}
					if($type eq "##attribute-ontology"){
						$headerInfo{$type}=$data[0]; #1 element
					}
					if($type eq "##species"){
						$headerInfo{$type}=$data[0]; #1 element
					}
					if($type eq "##genome-build"){
						$headerInfo{$type}=@data; #2 elements
					}
		        }
	        }
	    }
	    close($fh);
	}

    return \%headerInfo;
}

# @Purpose: Read a file from URL
# @input: 2 =>  String URL, String target (Target is not mandatory)
# @output: none
sub fetcher_JD {
	my ($url, $target) = @_;
    my $ua = LWP::UserAgent->new;
    $ua->timeout(10);
    $ua->env_proxy;

    my $response = $ua->get($url);
    if ($response->is_success) {
    	if($target){
	     	open my $OUT, '>', $target or die "File error: $! :: $?";
	        print $OUT $response->decoded_content;  # or whatever
	    }
	    else{
	    	my $string = $response->decoded_content;
	    	return $string ;
	    }
    }
    else {
        die $response->status_line;
    }
}

# @Purpose: retrieve the feature_ontology
# @input: 3 =>  String file, Hash, Int 
# @output: 1 => Object Ontology
sub _handle_ontology{
	my ($gff3headerInfo, $verbose, $quiet) = @_ ;

	my $ontology_obj=undef;
	my $internalO=1;
		
		if(exists_keys($gff3headerInfo, ("##feature-ontology"))){
			
			print "feature-ontology URI defined within the file: ".$gff3headerInfo->{'##feature-ontology'}."\n" if $verbose;		
			#retrieve the data from URI and save it in a string
			my $stringFILE=undef;
			try{
				$stringFILE = fetcher_JD($gff3headerInfo->{"##feature-ontology"});
			}
			catch{ 
				print "The URI provided (".$gff3headerInfo->{'##feature-ontology'}.") doesn't work.\n" if $verbose;
				print "error: $_\n" if ( $verbose >= 1);
			};

			if($stringFILE){
				#create a filehandler from a string
				open( my $fh_uriOnto, '<', \$stringFILE) || die "Cannot read the string: $! :: $?";
				
				#parse the ontology saved
		 		my $parser = undef;
		 		try{
		 			$parser = Bio::OntologyIO->new(-format => "obo",
	                                           -fh => $fh_uriOnto);
		 			$ontology_obj = $parser->parse();
		 			close $fh_uriOnto;
		 		}
		 		catch{
		 			print "The URI provided doesn't point to obo ontology format data.\n" if $verbose;
		 			print "error: $_\n" if ( $verbose >= 1);
		 			$parser = undef;
		 		};
		
				if($parser){ #We got ontology at the URI location, no need to use the internal one
					$internalO=undef;
					print "feature-ontology parsed correctly\n" if $verbose;
				}
			}
		}

	if($internalO){ #No URI provided for the feature-ontology(file case), or doesn't exist (hash / table case) let's use the interal one
	
		try{
			my $full_path = `perldoc -lm BILS::Handler::GXFhandler`;
			my $index = index($full_path, "BILS/");
			$index+=5; #To not chrinck BILS/ part of the path
			my $path_begin =  substr $full_path, 0, $index;
			my $correct_path = $path_begin."Ontology/SOFA";
			opendir (DIR, $correct_path) or die $!;
			my @list_file;
	
			# list all the sofa file available
			while (my $file = readdir(DIR)) {
				next if($file eq "." or  $file eq "..");
			 	push(@list_file, $file);
		    	}
	
		    	#get the most recent file
		    	my @sorted_list = sort { $a cmp $b } @list_file;
		    	my $recent_file = pop @sorted_list;
		    	my $sofa_file_path = $correct_path."/".$recent_file;
		    	print "We will use the most recent SOFA feature-ontology we have localy: $recent_file\n" if $verbose;

			#parse the ontology
			my $parser = Bio::OntologyIO->new(-format => "obo",
	                                      -file => $sofa_file_path);
			$ontology_obj = $parser->parse();
		}
		catch{
			print "error: $_\n" if ( $verbose >= 1);
 			print "Let's continue without feature-ontology information.\n" unless $quiet;
		};
	}
	return $ontology_obj;
}

1;
