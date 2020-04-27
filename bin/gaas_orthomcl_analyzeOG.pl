#!/usr/bin/env perl

use strict;
use warnings;
use Scalar::Util qw(looks_like_number);
#use Data::Dumper;
use Bio::DB::Taxonomy;
use Bio::TreeIO;
use Bio::Tree::Tree;
use Bio::Tree::TreeFunctionsI;
use Getopt::Long;
use IO::File;
use Pod::Usage;
use GAAS::GAAS;

my $header = get_gaas_header();

#VARIABLE DECLARATION
my $orthoMCL_file;
my $opt_help;
my $opt_output=undef;
my $outFile="output";
my $opt_tree;
my $nbProt=0;
my $speciesTreeString;
my $species_opt; my $focusThisTaxid="";
my $taxid_opt; my @TAXID_LIST;

my $message = "command line: orthomcl_analyzeOG.pl @ARGV\n\n";

# OPTION MANAGMENT
if ( !GetOptions( 'cog|kog|og=s'    => \$orthoMCL_file,
		          'o|out|output=s'  => \$opt_output,
                  't|tree=s'        => \$opt_tree,
                  'taxid=s'         => \$taxid_opt,
                  's|species=s'     => \$species_opt,
		          'h|help!'         => \$opt_help ) )

{
    pod2usage( { -message => "$header".'Failed to parse command line',
                 -verbose => 1,
                 -exitval => 1 } );
}

# Print Help and exit
if ($opt_help) {
    pod2usage( { -verbose => 99,
                 -exitval => 0,
                 -message => "$header\n" } );
}

if ( ! (defined($orthoMCL_file) ) ){
    pod2usage( {
           -message => "$header\nAt least 1 parameter is mandatory:\nInput reference KOG/COG file (--cog).\n\n",
           -verbose => 0,
           -exitval => 1 } );
}

my $outReport;
if ( $opt_output ){
    $outFile=$opt_output;
    $outFile=~ s/.gff//g;
    open($outReport, '>', $outFile."_report.txt") or die "Could not open file '$outFile' $!";
}
# print remind of options used
$message="Script launch ".localtime()."\n$message";
print $message; if ( $opt_output ){print $outReport $message;}
############################################

#### create connection to NCBI => useful to create tree / get taxid from species name
my $db = Bio::DB::Taxonomy->new(-source => 'entrez');
my %taxonList;

#### Manage others otpion:
my $focusThisSpecies;
if ( defined($species_opt) ){

    $focusThisTaxid="$species_opt";
    # get taxid
    if (! looks_like_number($focusThisTaxid)){ # try to retrive taxid from name
        $focusThisSpecies =  $focusThisTaxid;
        my @results=$db->get_taxonids($focusThisTaxid);
        $focusThisTaxid=$results[0];
    }# get scientific name
    else{ my $taxon  = get_taxon_efficiently($focusThisTaxid);
         $focusThisSpecies = $taxon->scientific_name; }

    $message = "You decided to focuse on $focusThisSpecies (taxid: $focusThisTaxid)\n";
    print $message; if ( $opt_output ){print $outReport $message;}
}

#### Manage output ######
my $outTree;   my $outTreeName=$outFile.'_gene_flux.nhx';
my $outSpTree; my $outSpTreeName=$outFile.'_species_tree.nhx';

if($opt_output){
    $outTree = Bio::TreeIO->new(
    -format => 'nhx',
    -file => '>'.$outTreeName,
    );

    $outSpTree = Bio::TreeIO->new(
    -format => 'nhx',
    -file => '>'.$outSpTreeName,
    );
}

my $screenDisplayTree = Bio::TreeIO->new(
    -format => 'nhx',
    -fh => \*STDOUT,
    );

##############
#### MAIN ####
##############

#### FIRST: PARSE GOs AND DEFINE TAXID LIST
$message = "Step: Parse GO File (Duplicate species by OG removed)\n";
print $message; if ( $opt_output ){print $outReport $message;}

my ($hashOG, $taxidListFromOG) = readOGfile($orthoMCL_file, $focusThisTaxid);

$message = "Parsing finished.\n";
print $message; if ( $opt_output ){print $outReport $message;}


######### get taxid from OG file if none given ###########
if ( ! (defined($taxid_opt) ) ){
    my $nbTaxid= @$taxidListFromOG;

    $message = "No species taxid list has been defined, we will use the $nbTaxid species present among the parsed OGs.\n";
    print $message; if ( $opt_output ){print $outReport $message;}

    @TAXID_LIST=@$taxidListFromOG;
}
else{ @TAXID_LIST=split(/[:,_\-\s\/]+/,$taxid_opt); }

$message = "\nList of taxid: \n";
print $message; if ( $opt_output ){print $outReport $message;}

foreach my $taxid (@TAXID_LIST){
    my $taxon  = get_taxon_efficiently($taxid);
    my $sci_name = $taxon->scientific_name;

    $message = "taxid $taxid corresponding to species $sci_name\n";
    print $message;     if ( $opt_output ){print $outReport $message;}
}

$message = "\n";
print $message;
if ( $opt_output ){print $outReport $message;}

######### SPECIES TREE MANAGEMENT ###########
my $speciesTreeReady;

$message = "Step: Species tree creation\n";
print $message; if ( $opt_output ){print $outReport $message;}

## Tree exist
if( defined ($opt_tree)){
    #Get Tree in File
    open(FIC,$opt_tree) or die "Couldn't open the file $opt_tree\n";
    #while( my $line = <FIC> ) {
    #    $line =~ s/\n//g;
    #    $speciesTreeString.=$line;
    #}
    my $test = Bio::TreeIO->new(-file => $opt_tree);
    $speciesTreeReady=$test->next_tree;
}
## Tree must be created
else{
    my $treeFromNCBI;
    my @species_names;
    foreach my $taxid (@TAXID_LIST){
        my $taxon = $db->get_taxon(-taxonid => "$taxid");
        my $spName=$taxon->scientific_name;
#       print "$taxid = $spName\n";
        push(@species_names, $spName);
    }
    $speciesTreeReady = $db->get_tree(@species_names);
#    print $speciesTreeReady, "\n";
    ## Clean Tree
    $speciesTreeReady->contract_linear_paths();
}

### Print Tree ###
$message = "This is the tree we will used:\n";
print $message; if ( $opt_output ){print $outReport $message;}

print $screenDisplayTree->write_tree($speciesTreeReady), "\n";
if ( $opt_output ){
    $outSpTree->write_tree($speciesTreeReady);

    #print species tree within the report file
    open F, "<$outSpTreeName" or die "Could not open file '$outSpTreeName' $!";;
        while (<F>) {
            print $outReport $_;
        }
    close F;
}

# create hash taxid => nodes (Leaves)
$message = "\n\nClean Taxid List according to tree provide:\n";
print $message; if ( $opt_output ){print $outReport $message;}

if (defined($opt_tree)){#Clean hashAbbTaxid to remove Taxid not existing in species Tree. Only useful when tree is no performed using taxid but comes as external file.
    my @copy_TAXID_LIST=@TAXID_LIST;
    foreach my $taxid (@copy_TAXID_LIST){
        my @nodes = $speciesTreeReady->find_node(-id => $taxid);
        if ($#nodes == -1){

            $message = "Species $taxid not present in Species Tree. We remove it from analyse.\n";
            print $message; if ( $opt_output ){print $outReport $message;}

            removeValueFromList($taxid, \@TAXID_LIST);
        }
        elsif ($#nodes >= 1){print "ERROR nb node not expected >1";exit;}
    }
}

## Create hast taxid/node whole Tree
my %hashTaxidNode;
foreach my $taxid (@TAXID_LIST){
    my @nodes = $speciesTreeReady->find_node(-id => $taxid);
    $hashTaxidNode{$taxid}=$nodes[0];
}

$message = "\nStep clean GO: keep only @TAXID_LIST\n";
print $message; if ( $opt_output ){print $outReport $message;}

my $hashOGFiltered = filterOGfileByTaxid($hashOG,\@TAXID_LIST);

$message = "\nStep: sort GO by species present \n";
print $message; if ( $opt_output ){print $outReport $message;}

my $hashOGsorted = sortOGforFlatDisplay($hashOGFiltered);
foreach my $key ( sort { $hashOGsorted->{$a} <=> $hashOGsorted->{$b}} keys %$hashOGsorted){
    $message = sizedPrint($key,100)."$hashOGsorted->{$key}\n";
    print $message;
    if ( $opt_output ){print $outReport $message;}
}

### Deduce appearance
$message = "\nStep: Deduce gene apparence (in ancestors or leaves)";
print $message; if ( $opt_output ){print $outReport $message;}

my $hashAppearance=deduceAppearance($hashOGFiltered, $speciesTreeReady, \%hashTaxidNode);
foreach my $taxID (keys %$hashAppearance){
    my $nbAppearance=$#{$hashAppearance->{$taxID}}+1;

    my $taxon  =  get_taxon_efficiently($taxID);
    my $sci_name = $taxon->scientific_name;

    $message = "$nbAppearance genes occured at $taxID node ($sci_name).\n";
    print $message; if ( $opt_output ){print $outReport $message;}

    my @node = $speciesTreeReady->find_node(-id => $taxID);
    $node[0]->add_tag_value('geneAppearance', $nbAppearance);
}

$message = "\nStep: Deduce gene losses\n";
print $message; if ( $opt_output ){print $outReport $message;}

### Deduce lost
my %LossIdByAppearance;
foreach my $taxID (keys %$hashAppearance){
    my @node = $speciesTreeReady->find_node(-id => $taxID);
    if( ! ($node[0]->is_Leaf )){ # avoid to study appearence on Leaf
        my %LossId;

        ## Create list of leaves
        my @original_ListLeavesId; my %original_hashNodeAllDescendant; my %original_hashTaxidAllDescendant;

        #print "\nAmong the $nbAppearance gene appeared at taxid $taxID ($sci_name) we have:\n";
        my @cladeNodes = $node[0]->get_all_Descendents(); ## I know all the descendent of Node of appearance
        foreach my $NodeFromDesc (@cladeNodes){ #get List leaves
            my $NodeFromDescTaxid = $NodeFromDesc->id();
            $original_hashNodeAllDescendant{$NodeFromDesc}=$NodeFromDescTaxid;
            $original_hashTaxidAllDescendant{$NodeFromDescTaxid}=$NodeFromDesc;
            if( $NodeFromDesc->is_Leaf ){
                push( @original_ListLeavesId, $NodeFromDescTaxid);
            }
        }
        ## Create the corresponding hash of leaves
        my %original_hashLeavesId = map { $_ => 1 } @original_ListLeavesId;
        ## Create the corresponding list of nodes
        my %original_hashLeavesNodes;
        foreach my $leafId (@original_ListLeavesId){
            my @node = $speciesTreeReady->find_node(-id => $leafId);
            $original_hashLeavesNodes{$leafId}=$node[0];
        }


        #### loop each list of present taxid
        my @ListsTaxidPresent = @{$hashAppearance->{$taxID}}; # all list of Leaves Present
        foreach my $oneList (@ListsTaxidPresent){
            my @ListTaxidPresent=@$oneList; # List of taxid that have the gene
            my %hashTaxidPresent = map { $_ => 1 } @ListTaxidPresent; # Hash of taxid that have the gene
#            print "ListTaxidPresent $#ListTaxidPresent @ListTaxidPresent  5555555555 sizeList $#original_ListLeavesId @original_ListLeavesId\n";


            ### Case No loss
            if($#ListTaxidPresent  == $#original_ListLeavesId){
                next;
            }

            ### Case only one loss
            elsif($#ListTaxidPresent  == $#original_ListLeavesId-1){
                 foreach my $taxid (@original_ListLeavesId){
                    if(! exists($hashTaxidPresent{$taxid})){ # avoid to study lost on Leaf
                        $LossId{$taxid}++;
                        last;
                    }
                }
            }

            ### Case several Lost or only one ancestral lost
            else{
                ## Create list of Taxid Absent
                my @ListTaxidAbsent;
                my %hashNodePresent; my %hashNodeAbsent;
                foreach my $taxid (@original_ListLeavesId){
                    my @node = $speciesTreeReady->find_node(-id => $taxid);
                    if(! exists($hashTaxidPresent{$taxid})){
                        push(@ListTaxidAbsent, $taxid);
                        $hashNodeAbsent{$node[0]}=$taxid;
                    }
                    else{$hashNodePresent{$node[0]}=$taxid;}
                }
                my %hashTaxidAbsent = map { $_ => 1 } @ListTaxidAbsent;

                ## Analyze
                my @listTaxidWithLoss;
                my %hashTaxidAdded;
                my @copy_listTaxidAbsent = @ListTaxidAbsent;
                while ($#copy_listTaxidAbsent > -1){
                    my $OneTaxid = shift @copy_listTaxidAbsent;
                    my $nodeOneTaxid = $original_hashTaxidAllDescendant{$OneTaxid};
#                    print "OneTaxid $OneTaxid  oneList List Absent: $@oneList = @ListTaxidAbsent\n";

                    ## Study ancestor
                    my $parentNode = $nodeOneTaxid->ancestor;
                    my @ancestorCladeNodes = $parentNode->get_all_Descendents();
                    my $presentGeneFound="no";
                    foreach my $anc_clNode (@ancestorCladeNodes){
                    #    my $clNodeTaxid = $clNode->id();
                    #    if(exists(clNodeTaxid))
                    #
                        if(exists($hashNodePresent{$anc_clNode})){
                            $presentGeneFound="yes";
                            $LossId{$OneTaxid}++;
#                           print "oui existe save child $OneTaxid\n";
                            last;
                        }
                    }
                    if($presentGeneFound eq "no"){
                        my $taxidFocused = $original_hashNodeAllDescendant{$parentNode};
                        if(! exists($hashTaxidAdded{$taxidFocused})) {
                            push ( @copy_listTaxidAbsent, $taxidFocused); # Push new ancestrak taxid to test if absent if before
                            $hashTaxidAdded{$taxidFocused}++;
#                            print "je push $taxidFocused \n";
                        }
                    }
                }
            }
        }
        # save results
        my $nbKey = keys %LossId;
        if ($nbKey >= 1){
            foreach my $key (keys %LossId){

                my $taxon  =  get_taxon_efficiently($key);
                my $sci_name = $taxon->scientific_name;
                #print "$LossId{$key} loss at $key ($sci_name)\n";
                $LossIdByAppearance{$taxID}{$key}=$LossId{$key};
            }
        }
        else{
            $LossIdByAppearance{$taxID}{'null'}++;
            #print "No loss\n";
        }
    }
}

## Merge All Loss A only one appearance
my %lossMergedByTaxid;
foreach my $keyID (keys %LossIdByAppearance){

    my $taxon  =  get_taxon_efficiently($keyID);
    my $sci_name = $taxon->scientific_name;
    my $nbAppearance=$#{$hashAppearance->{$keyID}}+1;

    my $message = "\nAmong the $nbAppearance genes appeared at taxid $keyID ($sci_name) we have:\n";
    print $message; if ( $opt_output ){ print $outReport $message; }

    if($LossIdByAppearance{$keyID}{'null'}){

        my $message = "No loss\n";
        print $message; if ( $opt_output ){ print $outReport $message; }
    }
    else{
        foreach my $keyID2 (keys %{$LossIdByAppearance{$keyID}}){
            my $value=$LossIdByAppearance{$keyID}{$keyID2};

            my $taxon  =  get_taxon_efficiently($keyID2);
            my $sci_name = $taxon->scientific_name;

            my $message = "$value loss at $keyID2 ($sci_name)\n";
            print $message; if ( $opt_output ){ print $outReport $message; }

            if(exists ($lossMergedByTaxid{$keyID2})){
                $lossMergedByTaxid{$keyID2}=$lossMergedByTaxid{$keyID2}+$value;
            }
            else{
                $lossMergedByTaxid{$keyID2}=$value;
            }
        }
    }
}

$message = "\nFinal losses Resume  (Total number of loss independent of their birth):\n";
print $message; if ( $opt_output ){ print $outReport $message; }

foreach my $keyID (keys %lossMergedByTaxid){

    my $taxon  =  get_taxon_efficiently($keyID);
    my $sci_name = $taxon->scientific_name;

    my $message =  "Gene lost in taxid $keyID ($sci_name) => $lossMergedByTaxid{$keyID}\n";
    print $message; if ( $opt_output ){ print $outReport $message; }

    my @node = $speciesTreeReady->find_node(-id => $keyID);
    $node[0]->add_tag_value('geneLoss', $lossMergedByTaxid{$keyID});
}

$message = "\n";
print $message; if ( $opt_output ){ print $outReport $message; }

print $screenDisplayTree->write_tree($speciesTreeReady), "\n";
if ( $opt_output ){
    print $outTree->write_tree($speciesTreeReady);

    #print tree within the report file
    open F, "<$outTreeName" or die "Could not open file '$outSpTreeName' $!";;
        while (<F>) {
            print $outReport $_;
        }
    close F;
}

my $finalMessage = "\nThe results given by this programm allow to have an overview of gene flux between the different species studied. The gene appearences and losses should be cautiously interpreted.".
"/!\\ The study is based on a Dollo like parsimomy. We allow genes to appear only one time. As you should be aware, annotations are often incomplete. Moreover, method of gene clustering have limitation to define precisely the orthologous groups.".
" Consequently, the gene described as lost here are only potential losses. Verification should be performed.\nEND\n";

print $finalMessage;
if ( $opt_output ){ print $outReport $finalMessage; }

################################################################## FUNCTIONS #####################################
sub deduceAppearance{
    my ($hashOGFiltered, $speciesTreeReady, $hashTaxidNode)=@_;

    my %hashAppearance;
    my $nbOGstudied=0;
    foreach my $key (keys %$hashOGFiltered){
        my @speciesList=@{$hashOGFiltered->{$key}};
        if ($#speciesList == 0){
                $nbOGstudied++;
                push (@{$hashAppearance{$speciesList[0]}}, [@speciesList]);
        }
        else{
            my @nodesList;
            foreach my $taxid (@speciesList){
                    push (@nodesList, $hashTaxidNode{$taxid});
            }
            my $nbNodes=scalar @nodesList;
    #       print "nbNodes $nbNodes\n";
            if ( $nbNodes <= 1){
    #           print "Not enough species kept for reconstruct lca. We skip this OG.\n";
                next;
            }
            my $lca=$speciesTreeReady->get_lca(-nodes => \@nodesList);
            my $idLCA=$lca->id();
            $nbOGstudied++;
     #       print "resu $nbOGstudied $idLCA\n";
            push (@{$hashAppearance{$idLCA}}, [@speciesList]);
        }
    }
    print "\nWe studied $nbOGstudied OGs\n";
    return (\%hashAppearance);
}

sub removeValueFromList{
    my ($val, $array)=@_;
    my $index = 0;
    $index++ until $array->[$index] eq $val;
    splice(@$array, $index, 1);
}

sub sortOGforFlatDisplay{
my ($hashOGref)=@_;
my %hashOGidSentence;
    foreach my $OGkey (keys %$hashOGref){
        my @ListSpeciesOG=@{$hashOGref->{$OGkey}};
        my @ListSpeciesOGSorted =  ( sort ({ $a <=> $b } @ListSpeciesOG));
        my $IDsentence; my $cpt=0;
        foreach my $key (@ListSpeciesOGSorted){
            if ($cpt == 0){$IDsentence.="$key";}
            else{$IDsentence.="_$key";}
            $cpt++;
        }
        $hashOGidSentence{$IDsentence}++;
    }
    return \%hashOGidSentence;
}

sub filterOGfileByTaxid{
my ($hashOGref, $ListTaxidToTest)=@_;
my $OGnotKept=0;
my $OGKept;
my %hashOGrefCleaned;
my %taxidAnalyzed;

    foreach my $OGkey (keys %$hashOGref){ # foreach OG group
        my @newOG;
        my @ListSpeciesOG=@{$hashOGref->{$OGkey}};
        foreach my $taxidOG (@ListSpeciesOG){ # foreach species of the OG
            foreach my $taxidToTest (@$ListTaxidToTest){
                if($taxidOG eq $taxidToTest){
                    push(@newOG, $taxidOG);
                    $taxidAnalyzed{$taxidOG}++;
                }
            }
        }
        if ($#newOG == -1 ){
            $OGnotKept++;
        }
        else{
            @{$hashOGrefCleaned{$OGkey}}=@newOG;
            $OGKept++;
        }
    }
    print "$OGnotKept OG removed\n"; print "$OGKept OG kept containing at least one of these species: @$ListTaxidToTest\n";

    ##check if Taxid in list given are not present in OG => Warning message
    foreach my $taxidToTest (@$ListTaxidToTest){
        if ( ! exists ($taxidAnalyzed{$taxidToTest})){ print "///!\\\\\\ WARNING MESSAGE: Taxid $taxidToTest given as input is not present in the OG analyzed."}
    }

    return \%hashOGrefCleaned;
}

sub readOGfile{
my ($OG_file, $TaxidThisSpecies)=@_;
my %allTaxidFromOG; my @allTaxid;
my %nbSpeciesByOG;
my %nbSpeciesByOGwithThisSpecies;
my %HashSpeciesByOG;

    open(FIC,$orthoMCL_file) or die "Couldn't open the file $orthoMCL_file\n";
    while( my $line = <FIC> ) {
      	if (! ($line =~ /^#/)){
            chomp($line) ;
            my @splitedLine = split(" ",$line);
            my %hashOGspecies; my %hashOGspeciesWithThisSpecies; my $OGcontainsThisSpecies="no";

            my $OGname=shift @splitedLine;
            $OGname=~s/://g ; # remove ":"
            foreach my $prot (@splitedLine ){
                my @splitedProt = split("\\|",$prot);
                $hashOGspecies{$splitedProt[0]}++; #hash of taxid
                $allTaxidFromOG{$splitedProt[0]}++;
                # check for species centered analysis
                if ($splitedProt[0] eq $TaxidThisSpecies){
                    $OGcontainsThisSpecies="yes";
                }
            }
            # species centered analysis
            if ($OGcontainsThisSpecies eq "yes"){
                foreach my $prot (@splitedLine ){
                    my @splitedProt = split("\\|",$prot);
                    $hashOGspeciesWithThisSpecies{$splitedProt[0]}++; #hash of taxid
                }
                my $nbSpecies= keys %hashOGspeciesWithThisSpecies; # count number of taxid
                $nbSpeciesByOGwithThisSpecies{$nbSpecies}++; #General infornation about number of OG containing a number of species
            }

            my $nbSpecies= keys %hashOGspecies; # count number of taxid
            $nbSpeciesByOG{$nbSpecies}++; #General infornation about number of OG containing a number of species
        #   print "There is $nbSpecies species in the group $OGname\n";
            foreach my $key (keys %hashOGspecies){
                push(@{$HashSpeciesByOG{$OGname}}, $key);
            }
        }
    }

    ### Manage list of all taxid present in OGs
    foreach my $key (keys %allTaxidFromOG){
        push (@allTaxid, $key)
    }

    my $nbOGanalyzed=keys %HashSpeciesByOG;
    print "$nbOGanalyzed OGs analized:";

    foreach my $nbKey (sort {$b <=> $a} keys %nbSpeciesByOG){
    	print "There is $nbSpeciesByOG{$nbKey} OG containing $nbKey species\n";
    }
    print "\n\n";
    if($TaxidThisSpecies ne ""){
        print "Analysis focusing on $TaxidThisSpecies\n";
        my $size = keys %nbSpeciesByOGwithThisSpecies;
        if ($size == 0){print "This taxid has not been retrieved among GOs parsed. \n";}
        else{
            foreach my $nbKey (sort {$b <=> $a} keys %nbSpeciesByOGwithThisSpecies){
                print "There is $nbSpeciesByOGwithThisSpecies{$nbKey} OG containing $nbKey species including taxid $TaxidThisSpecies\n";
            }
        }
    }
 return \%HashSpeciesByOG, \@allTaxid ;
}

sub get_taxon_efficiently{
    my ($taxid) = @_;

    my $taxon;
    if(! exists($taxonList{$taxid})){
        if($taxid !~ /^\d+$/){warn "taxid is expected to be a number. We got <$taxid>. Please fix your txt file to have proper taxid.";exit;}
        $taxon = $db->get_taxon(-taxonid => $taxid);
        $taxonList{$taxid} = $taxon;
    }else{
        $taxon = $taxonList{$taxid};
    }

}
__END__

=head1 NAME

analyzeOG.pl -
The script computes some statistics of a COG/KOG file from OrthoMCL output -
Statistics as : - number of OG by number of species
                - number of OG by number of species that includes a specifc species (if specified by -s option)
                - gene appearances
                - gene losses
                - Amount of gene at each node/leaf (Compute for ancestral nodes is Dallo like parsimony => on gene appear only one time. Moreover we do not consider potential HGTs).
                - etc...
*OG = Ortholog Group

Prerequisite: OrthoMCL output where taxid (http://www.ncbi.nlm.nih.gov/taxonomy/) has been used to label the sequences. As example:
OG00001: 10090|ENSMUSP0000001 9606|ENSP0000001

=> In that case The orthologous group "OG00001" contains 2 sequences, one comming from Mus musculus (taxid 10090) and the other from human (taxid 9606).

=head1 SYNOPSIS

    ./analyzeOG.pl --cog infile [ -s TaxidSpeciesWhichFocusOn -t tree -taxid taxidA_taxidB_taxidC --output outfile ]
    ./analyzeOG.pl --help

=head1 OPTIONS

=over 8

=item B<--cog>,  B<--og> or B<--kog>

Orthomcl file containg Ortholog groups (COG) from OrthoMCL.

=item B<--taxid>

Taxid list. If provided the analyse will use only these species. If a tree is also provided, the taxid will be filtered according to the tree to keep only taxid present in the tree.
If no taxid is provided, but a tree is, only species from the tree will be analyzed. If no tree and no taxid are provided, only taxid among OG will be use.

=item B<-t> or B<--tree>

Tree file in nhx format. If provided the analyse will focuse only on species present in the tree.
When no tree is provided, a species tree will be created on the fly using the NCBI taxonomy database online according to the species present among the OG.

=item B<-s>, B<--species> or B<-ref>

taxid or scientific name (use underscore instead of spaces). It allows to focus the analysis only on OG containg the species defined.

=back

=cut
