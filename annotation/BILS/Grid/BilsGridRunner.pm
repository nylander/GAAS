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
    ## find the path to the SLURM_perl_lib directory.
    foreach my $dir (@INC) {
        if (-d "$dir/SLURM_perl_lib") {
            push (@INC, "$dir/SLURM_perl_lib");
            last;
        }
    }
}

use Run_Bsub;
use Cwd;

####
sub run_on_grid {
    my ($args) = @_  ;

    my (@cmds, $scheduler, $queueu);
    if( ! defined($args->{cmds}) ) {print "No command provided.\n";exit;}    		 else{ @cmds = $args->{cmds}; }
    if( ! defined($args->{scheduler})) {$scheduler="slurm"; print "Default scheduler used: $scheduler\n"} 		 else{ $scheduler = $args->{scheduler}; }
    if( ! defined($args->{queueu})) {$queueu=undef;} 		 else{ $queueu = $args->{queueu}; }

    @cmds = shuffle @cmds;

    if($scheduler eq "slurm"){
      &Run_Slurm::set_queue("normal");
      &Run_Slurm::set_memory("4"); # 4 G of RAM
      #&Run_Bsub::set_mount_test(cwd()); # only run on nodes with a verified mount
      my @failed_cmds = &Run_Slurm::run(@cmds);
    }
    elsif($scheduler eq "lsf"){
      &Run_Bsub::set_queue("normal");
      &Run_Bsub::set_memory("4"); # 4 G of RAM
      #&Run_Bsub::set_mount_test(cwd()); # only run on nodes with a verified mount
      my @failed_cmds = &Run_Bsub::run(@cmds);
    }
    else{
        print "Scheduler $scheduler not implemented yet. <slurm> and <lsf> are the only possible choice currently\n";
    }

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
