use strict;
use warnings;

use LWP::Simple;
use File::Fetch;

	my $sofa_url="https://raw.githubusercontent.com/The-Sequence-Ontology/SO-Ontologies/master/subsets/SOFA.obo";
	my $ff = File::Fetch->new(uri => $sofa_url);
	#my $sofa_url="https://github.com/The-Sequence-Ontology/SO-Ontologies/blob/master/subsets/SOFA.obo";
	my $sofa_file_path = $path."sofa.obo";
	if(! -f $sofa_file_path){
		print "SOFA file missing, we download it\n";

		eval {
 			getstore($sofa_url, $sofa_file_path);
 			}; if ($@) {
 				print "We got an error: $@\n";
 				die;
 			}
 		my $file = $ff->fetch(to => $sofa_file_path);
 		print Dumper($file);
	}
	exit;
