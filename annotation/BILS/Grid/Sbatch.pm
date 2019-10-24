package Sbatch;

use strict;
use warnings;
use File::Basename;
#use IPC::Cmd qw[can_run run];
use Carp;
use Moose;
use BILS::Grid::GridRunner;

extends 'GridRunner';

has scheduler => ('is' => 'rw', isa => 'Str', default => 'Slurm');

# The BUILD method is called after an object is created.
sub BUILD {
     my $answer = system("squeue 1>/dev/null 2>&1");
     if ($? != 0) {
         print "Slurm does not seem to be installed!\n";  exit;
    }
}


####
sub _submit_job {
    my $self = shift;
    my $num_cmds_launched = shift;

    my $num_cmds = $self->{num_cmds};
    my $cmds_list_aref = $self->{cmds_list};

    my $log_dir = $self->{log_dir};
    my $retvals_dir = $self->{retvals_dir};
    my $cmds_dir = $self->{cmds_dir};
    my $monitor_dir = $self->{monitor_dir};

    my $orig_num_cmds_launched = $num_cmds_launched;


    # ---------------- create the job (write into a shell script) ---------------
    my $shell_script = "$cmds_dir/J$$.S${num_cmds_launched}.sh";
    open (my $fh, ">$shell_script") or die $!;
    print $fh "#!/bin/sh\n\n";
    $self->_write_minimal_environment($fh);

    my $monitor_started = "$monitor_dir/$num_cmds_launched.started";
    my $monitor_finished = "$monitor_dir/$num_cmds_launched.finished";

    my @cmd_indices_prepped;

    my $next_cmd_index = $num_cmds_launched; #always one less than the current index
    my $cmd_string = $cmds_list_aref->[ $next_cmd_index ];

    push (@cmd_indices_prepped, $next_cmd_index);

		my $retval_bin = int($next_cmd_index / $self->retval_bin_size());

		my $retval_subdir = "$retvals_dir/$retval_bin";
		unless (-d $retval_subdir) {
			mkdir $retval_subdir or die "Error, cannot mkdir $retval_subdir";
		}

    print $fh "## Command index $next_cmd_index\n"
        . "touch $monitor_started\n"
        . "$cmd_string\n"
        . 'echo $? >> ' . "$retval_subdir/entry_$next_cmd_index.ret\n\n";

    print $fh "\n"
        . "rm -f $monitor_started\n"
        . "touch $monitor_finished\n"
        . "\n"
        . "exit 0\n\n";


    close $fh;
    chmod (0775, $shell_script);


    # ---------------- create the bsub command ---------------
    my $cmd = undef;
    if($self->{queue}){
      my $queue = $self->{queue};
      $cmd = "sbatch -p $queue -e $shell_script.stderr -o $shell_script.stdout ";
    }
    else{
	$cmd = "sbatch -e $shell_script.stderr -o $shell_script.stdout ";
    }
    if (my $memory = $self->{memory}) {
  	$cmd .= " --mem=".$memory."gb ";
    }
    $cmd .= " $shell_script 2>&1 ";

    # ---------------- run the bsub job ---------------
    print "Submitting: $shell_script with sbatch\n" if $self->verbose;
    my $job_id_text = `$cmd`;
    $num_cmds_launched++;

    # ---------------- check status ---------------
    my $ret = $?;
    if ($ret) {
        print STDERR "SBATCH failed to accept job: $cmd\n (ret $ret)\n";

        unlink $shell_script; # cleanup, try again later

        sleep(2*60); # sleep 2 minutes for now.  Give the system time to recuperate if a problem exists
        return ($orig_num_cmds_launched);

    }
    else {

        $shell_script = basename($shell_script);
        open (my $logdir_jobsfh, ">>$log_dir/job_ids.txt") or die "Error, cannot open file $log_dir/job_ids.txt";
        ## get the job ID and log it:
        if ($job_id_text =~ /Submitted batch job (\d+)/) {
            my $job_id = $1;

            print $logdir_jobsfh "$job_id\t$shell_script\n";

            $self->{nodes_in_progress}->{$monitor_finished} = $job_id;
            $self->{job_id_to_cmd_indices}->{$job_id} = \@cmd_indices_prepped;
            $self->{job_id_to_submission_time}->{$job_id} = time();
        }
        else {

            die "Fatal error, couldn't extract Job ID from submission text: $job_id_text";

        }
        close $logdir_jobsfh;

        # sleep($WAITTIME); # wait just a short while to give the system a few seconds to act on the submitted jobs.
        return ($num_cmds_launched);
    }

}

####
sub _job_running_or_pending_on_grid {
    my $self = shift;
    my ($job_id) = @_;

    if (time() - $self->{job_id_to_submission_time}->{$job_id} < $self->resort_to_polling_time) {
        return("TOO_SOON");
    }


    # print STDERR "Polling grid to check status of job: $job_id\n";

    my $response = `squeue $job_id`;
    #print STDERR "Response:\n$response\n";

    foreach my $line (split(/\n/, $response)) {
        my @x = split(/\s+/, $line);

        if ($x[0] eq $job_id) {
            my $state = $x[2];
            if ($state eq "DONE" || $state eq "EXIT") {
                return(0);
            }
            else {
                $self->{job_id_to_submission_time}->{$job_id} = time();
                return($state);
            }
        }
    }

    print STDERR "-no record of job_id $job_id, setting as state unknown\n";
    return undef; # no status info
}

#no Moose;
#__PACKAGE__->meta->make_immutable;
;
1
