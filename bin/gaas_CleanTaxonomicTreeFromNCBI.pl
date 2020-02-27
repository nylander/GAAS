#!/usr/bin/perl

use strict;
use warnings;
use Bio::DB::Taxonomy;
use Bio::TreeIO;
use Bio::Tree::NodeNHX;
use Getopt::Long;
use IO::File;
use Pod::Usage;

#VERIABLE DECLARATION

my $opt_tree;
my $opt_help;
my $opt_output;
my $nbProt=0;

# OPTION MANAGMENT
if ( !GetOptions( 't=s' => \$opt_tree,
				  'o|output=s'      => \$opt_output,
				  'h|help!'         => \$opt_help ) )

{
    pod2usage( { -message => 'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

if ($opt_help) {
    pod2usage( { -verbose => 2,
                 -exitval => 0 } );
}
 
if ( ! ( (defined($opt_tree)) ) ){
    pod2usage( {
           -message => "\nAt least 1 parameter is mandatory:\nInput reference file (--t)\n\n",
           -verbose => 0,
           -exitval => 2 } );
}

#my $out = new Bio::TreeIO->new(-file => '>'.$opt_output.'.svg',
#                          -format => 'svggraph');


#Get Tree in File
my $treeString;
open(FIC,$opt_tree) or die "Couldn't open the file $opt_tree\n";
while( my $line = <FIC> ) {
	$line =~ s/\n//g;
	$treeString.=$line;
}

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
			$comaBetween="yes";print "We found >>>, \n";
		}

		if($char eq ")"){
			print "We found >>>) \n";
			$nbClosedP++;
			if($nbClosedP == 1){$positionLastP=$endPosition;}
			
			if($nbClosedP > 1){
				
				#define position of two last closed parenthesis
				if($nbClosedP == 2){$positionBeforeLastP=$endPosition;}
				if($nbClosedP > 2){$positionLastP=$positionBeforeLastP; $positionBeforeLastP=$endPosition;}
				
				### Unitary Test (XXXX)blabla
				print "Unitary Test comaBetween $comaBetween\n";
				my $currentinfo=substr $currentTString, $positionLastP;
				my $positionUniq=checkUniqParenthesis($positionLastP-1, $currentTString); #return "("
				if ($positionUniq != -1) {
					print "We have to remove \n";
					my $positionEndCut=getPositionEndCut($positionLastP+1,$currentTString);
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
					$comaBetween="no";last;
				}
				
				### CHECK opening parenthesis // get position
				if ($comaBetween eq "no"){
					print "Second Test\n";
					my $posiotionFirstOP=checkOpeningParenthesis($positionLastP, $currentTString);
					my $posiotionAfterFirstOP=checkOpeningParenthesis($positionBeforeLastP, $currentTString);
					
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
						print "C=>newStringPart1 $newStringPart1\n";
						my $positionEndCut=getPositionEndCut($positionLastP,$newTString);
						$newStringPart2=substr $newTString, $positionEndCut;
						$newTString="$newStringPart1$newStringPart2";
						print "C=>newStringPart2 $newStringPart2\n";
						$sizeAfter=length($newTString);
#						print "REMOVING closing parenthesis =>new Size = $sizeAfter \n\n";
						$comaBetween="no";last;
					}		
				}
			}
			$comaBetween="no";
		}
		$endPosition--;
	}
	$sizeTreeBefore=$sizeCurrentTree; # print "sizeTreeBefore $sizeTreeBefore\n";
	$sizeCurrentTree=length($newTString); # print "sizeCurrentTree $sizeCurrentTree
					print"\nEND OF ROUND !!$newTString\n\n"; 
	$currentTString=$newTString;
}

### final case / last round:
### Unitary Test (XXXX)blabla
my $positionUniq=checkUniqParenthesis($positionBeforeLastP-1, $currentTString); #return "("
if ($positionUniq != -1) {
	my $positionEndCut=getPositionEndCut($positionBeforeLastP+1,$currentTString);
	my $newStringPartA=substr $currentTString, 0, $positionUniq;
	my $newStringPartB=substr $currentTString, $positionUniq+1;
	my $currentTStringMinusPO="$newStringPartA$newStringPartB";
	$newStringPart1=substr $currentTStringMinusPO, 0, $positionBeforeLastP-1;
	$newStringPart2=substr $currentTStringMinusPO, $positionEndCut-1;
	$newTString="$newStringPart1$newStringPart2";
}

print "\nfinalTree= $newTString\n";

#open(my $fh, '>', "treeCleanResult.txt") or die "Could not open file 'treeCleanResult' $!";
#print $fh "$finalTree\n";
#close $fh;



##### METHODS #######

sub getPositionEndCut{
	my ($endPosition,$treeString)=@_;
#	print "getPositionEndCut $endPosition\n";
	while($endPosition < length($treeString)){
		my $char = substr $treeString, $endPosition, 1;
		print "ENDCUT char = $char \n";
		if($char eq ',' || $char eq ')' || $char eq ';'){
			return $endPosition;
		}
	$endPosition++;
	}
}

sub checkOpeningParenthesis{
	my ($positionClosedP, $treeString)=@_;
	my $endPosition=$positionClosedP;
	my $positionOP=0;
	my $nbPopened=0; my $nbPclosed=0;
	while($endPosition >= 0){
		my $char = substr $treeString, $endPosition, 1;
#		print "char = $char $endPosition <= $nbPopened $nbPclosed\n";
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

sub checkUniqParenthesis{
	my ($endPosition, $treeString)=@_;
	print "checkUniqParenthesis $endPosition\n";
	my $nbPclosed=0;
	while($endPosition >= 0){
		my $char = substr $treeString, $endPosition, 1;
		print "charUniq = $char $endPosition\n";
		if( ( $char eq ')') || ($char eq ",") ){ return -1;}
		if($char eq '('){
			return $endPosition;		
		}
		$endPosition--;
	}
}

