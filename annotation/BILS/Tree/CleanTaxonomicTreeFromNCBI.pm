#!/usr/bin/perl -w

package BILS::Tree::CleanTaxonomicTreeFromNCBI;

use Data::Dumper;


use base qw/Exporter/;

our @EXPORT = qw/ clean_taxonomic_tree_from_ncbi /;


our @EXPORT_OK = qw/ some_other /;

=head1 SYNOPSIS


=head1 DESCRIPTION

	A library to clean the taxonomic tree coming from NCBI. Indeed there is several ancestor in some unitary internal branches
	Inherits from BILS::Tree
	
=cut	



sub clean_taxonomic_tree_from_ncbi{
	my ($treeString)=@_;

	# Analyse the TREE string
	my $newTString=$treeString;
	my $currentTString=$treeString;
	my $sizeTreeBefore=0;
	my $sizeCurrentTree=length($treeString);
	my $newStringPart1; my $newStringPart2;
	my $positionLastP; my $positionBeforeLastP;

	while( ! ($sizeTreeBefore == $sizeCurrentTree)){

		#parcours de l'abre
		my $endPosition=$sizeCurrentTree-1;
		my $nbClosedP=0; my $comaBetween="no";	

		while($endPosition != 0){
			my $char = substr $currentTString, $endPosition, 1;
			
			if (($char eq ",") && ($nbClosedP>=1)){
				$comaBetween="yes";
			}

			if($char eq ")"){
				$nbClosedP++;
				if($nbClosedP == 1){$positionLastP=$endPosition;}
				
				if($nbClosedP > 1){
					
					#define position of two last closed parenthesis
					if($nbClosedP == 2){$positionBeforeLastP=$endPosition;}
					if($nbClosedP > 2){$positionLastP=$positionBeforeLastP; $positionBeforeLastP=$endPosition;}
					
					### Unitary Test (XXXX)blabla
					my $positionUniq=__checkUniqParenthesis__($positionLastP-1, $currentTString); #return "("
					if ($positionUniq != -1) {
						my $positionEndCut=__getPositionEndCut__($positionLastP+1,$currentTString);
						### remove opening parenthesis part
						my $newStringPartA=substr $currentTString, 0, $positionUniq;
						my $newStringPartB=substr $currentTString, $positionUniq+1;
						my $currentTStringMinusPO="$newStringPartA$newStringPartB";
						### remove closing parenthesis part
						$newStringPart1=substr $currentTStringMinusPO, 0, $positionLastP-1;
	#					print "A=>newStringPart1 $newStringPart1\n";
						$newStringPart2=substr $currentTStringMinusPO, $positionEndCut-1;
	#					print "A=>newStringPart2 $newStringPart2\n";
						$newTString="$newStringPart1$newStringPart2";
						last;
					}
					
					### CHECK opening parenthesis // get position
					if ($comaBetween eq "no"){
						my $posiotionFirstOP=__checkOpeningParenthesis__($positionLastP, $currentTString);
						my $posiotionAfterFirstOP=__checkOpeningParenthesis__($positionBeforeLastP, $currentTString);
						
						### IF closeing parenthesis are in consecutive position
						if(($posiotionAfterFirstOP) == ($posiotionFirstOP+1)){ 						
							#### remove opening parenthesis part // comsequences shift position to minus 1 in the newTString
	#						print "CONSECUTIVE POSITION FOR OPENING PARENTHESIS =>  \n";
							$newStringPart1=substr $currentTString, 0, $posiotionFirstOP;
							#print "B=>newStringPart1 $newStringPart1\n";
							my $posiotionGetEnd=$posiotionFirstOP+1;
							$newStringPart2=substr $currentTString, $posiotionGetEnd;
							#print "B=>newStringPart2 $newStringPart2\n";					
							my $sizeBeforenewString=length($newTString);
							$newTString="$newStringPart1$newStringPart2"; # create the new string
							my $sizeAfter=length($newTString);

							#### remove closed parenthesis part of the new string
							$newStringPart1=substr $newTString, 0, $positionLastP-1;
	#						print "C=>newStringPart1 $newStringPart1\n";
							my $positionEndCut=__getPositionEndCut__($positionLastP,$newTString);
							$newStringPart2=substr $newTString, $positionEndCut;
							$newTString="$newStringPart1$newStringPart2";
	#						print "C=>newStringPart2 $newStringPart2\n";
							$sizeAfter=length($newTString);
	#						print "REMOVING closing parenthesis =>new Size = $sizeAfter \n\n";
							last;
						}		
					}
				}
				$comaBetween="no";
			}
			$endPosition--;
		}
		$sizeTreeBefore=$sizeCurrentTree; # print "sizeTreeBefore $sizeTreeBefore\n";
		$sizeCurrentTree=length($newTString); # print "sizeCurrentTree $sizeCurrentTree\nEND OF ROUND !!$newTString\n\n"; 
		$currentTString=$newTString;
	}

	### final case / last round:
	### Unitary Test (XXXX)blabla
	my $positionUniq=__checkUniqParenthesis__($positionBeforeLastP-1, $currentTString); #return "("
	if ($positionUniq != -1) {
		my $positionEndCut=__getPositionEndCut__($positionBeforeLastP+1,$currentTString);
		my $newStringPartA=substr $currentTString, 0, $positionUniq;
		my $newStringPartB=substr $currentTString, $positionUniq+1;
		my $currentTStringMinusPO="$newStringPartA$newStringPartB";
		$newStringPart1=substr $currentTStringMinusPO, 0, $positionBeforeLastP-1;
		$newStringPart2=substr $currentTStringMinusPO, $positionEndCut-1;
		$newTString="$newStringPart1$newStringPart2";
	}
	return $newTString;
}


##### METHODS #######

sub __getPositionEndCut__{
	my ($endPosition,$treeString)=@_;
#	print "getPositionEndCut $endPosition\n";
	while($endPosition < length($treeString)){
		my $char = substr $treeString, $endPosition, 1;
		if($char eq ',' || $char eq ')' || $char eq ';'){
			return $endPosition;
		}
	$endPosition++;
	}
}

sub __checkOpeningParenthesis__{
	my ($positionClosedP, $treeString)=@_;
	my $endPosition=$positionClosedP;
	my $positionOP=0;
	my $nbPopened=0; my $nbPclosed=0;
	while($endPosition >= 0){
		my $char = substr $treeString, $endPosition, 1;
		if($char eq '('){$nbPopened++;}
		if($char eq ')'){$nbPclosed++;}
		if($nbPopened eq $nbPclosed){
			$positionOP=$endPosition;
			return $positionOP;
		}
		$endPosition--;
	}
	print "No symetric open parenthesis !\n";exit;
}

sub __checkUniqParenthesis__{
	my ($endPosition, $treeString)=@_;
#	print "checkUniqParenthesis $endPosition\n";
	my $nbPclosed=0;
	while($endPosition >= 0){
		my $char = substr $treeString, $endPosition, 1;
		if( ( $char eq ')') || ($char eq ",") ){ return -1;}
		if($char eq '('){
			return $endPosition;		
		}
		$endPosition--;
	}
}

1;