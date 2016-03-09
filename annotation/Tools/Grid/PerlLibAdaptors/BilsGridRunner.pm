package BilsGridRunner;

use strict;
use warnings;
use Carp;
use base qw (GridRunner);  ## really just for aesthetics
use List::Util qw (shuffle);
use FindBin;

BEGIN {

    ## find the path to the LSF_perl_lib directory.
    foreach my $dir (@INC) {
        if (-d "$dir/LSF_perl_lib") {
            push (@INC, "$dir/LSF_perl_lib");
            last;
        }
    }
    
}

use Run_Bsub;
use Cwd;

####
sub run_on_grid {
    my @cmds = @_;

    @cmds = shuffle @cmds;

    &Run_Bsub::set_queue("normal"); 
    &Run_Bsub::set_memory("4"); # 4 G of RAM
    #&Run_Bsub::set_mount_test(cwd()); # only run on nodes with a verified mount
    
    my @failed_cmds = &Run_Bsub::run(@cmds);

    if (@failed_cmds) {
        
        my $num_failed_cmds = scalar(@failed_cmds);
        
        print STDERR "$num_failed_cmds commands failed during grid computing.\n";

        return(&run_on_grid(@failed_cmds));
        
    }
    else {
        print "All commands completed successfully on the computing grid.\n";
        return(0);
    }
}

####

    
1;

    

