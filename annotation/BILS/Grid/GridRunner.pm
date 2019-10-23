package GridRunner;

use strict;
use warnings;
use Cwd;
use Carp;
use vars qw ($log_dir $keep_logs_on_failure);
use Moose;



# define attributes
has waitime => ('is' => 'rw', isa => 'Int', default => 15 );
has retval_bin_size => ('is' => 'rw', isa => 'Int', default => 1000 );
has resort_to_polling_time => ('is' => 'rw', isa => 'Int', default => 15*60 );# 15 minutes
has job_id_to_prevtime => ('is' => 'rw', isa => 'HashRef');
has max_job_at_a_time => ('is' => 'rw', isa => 'Int', default => 500 );
has queue => ('is' => 'rw', isa => 'Str');
has memory => ('is' => 'rw', isa => 'Int', default => 4);
has cmds_list => ('is' => 'rw', isa => 'ArrayRef', required => 1);
has num_cmds => ('is' => 'rw', isa => 'Int');
has log_dir => ('is' => 'rw', isa => 'Str');
has cmds_dir => ('is' => 'rw', isa => 'Int', default => 4);
has retvals_dir => ('is' => 'rw', isa => 'Int', default => 4);
has monitor_dir => ('is' => 'rw', isa => 'Int', default => 4);
has retvalues => ('is' => 'rw', isa => 'Int', default => 0); #later set to array ref with retvalues for each command.
has nodes_in_progress => ('is' => 'rw', isa => 'HashRef' );
has job_id_to_cmd_indices => ('is' => 'rw', isa => 'Hash' ); # job_id => [1,2,3,4,...]  so we know which cmds correspond to each grid job identifier.
has job_id_to_submission_time => ('is' => 'rw', isa => 'Hash' );
has mount_test => ('is' => 'rw');
has verbose => ('is' => 'rw', isa => 'Bool', default => 0);



# The BUILD method is called after an object is created.
# Here it is used to set all different folders used to store logging
sub BUILD {
    my $self = shift;

# -------- set up num_cmds -----------

    $self->{num_cmds} = scalar @{$self->cmds_list};

# --------------------------------------------
# -------- Set up log environment ------------
# --------------------------------------------

    # -------- create a log directory ------------
    my $log_dir;
    unless ($self->log_dir) {

        # create logging area
        $log_dir = cwd;
        unless (-d $log_dir) {
            mkdir ($log_dir);
        }
    }

    # -------- create a log sub_directory ------------
    my $hostname = `hostname`;
	  chomp($hostname);
    $log_dir .= "/job.J$$.$hostname.$$." . time();
    unless (-d $log_dir) {
        mkdir $log_dir or die "Error, cannot mkdir $log_dir";
    }

    # --------- write commands listing ------------
    open (my $fh, ">$log_dir/cmds_list.txt") or die $!;
    my $index = 0;
    foreach my $cmd ($self->cmds_list) {
        print $fh "index($index)\t$cmd\n";
        $index++;
    }
    close $fh;


    # ------------ finish logdir setup and object creation. ------------
    my $cmds_dir = "$log_dir/cmds";
    my $retvals_dir = "$log_dir/retvals";
    my $monitor_dir = "$log_dir/monitor";
    foreach my $dir ($cmds_dir, $retvals_dir, $monitor_dir) {
        mkdir $dir or die "Error, cannot mkdir $dir";
    }

    $self->{log_dir} = $log_dir,
    $self->{cmds_dir} = $cmds_dir,
    $self->{retvals_dir} = $retvals_dir,
    $self->{monitor_dir} =  $monitor_dir,
}


sub _get_num_nodes_used {
    my $self = shift;
    my $num_nodes_used = scalar (keys %{$self->{nodes_in_progress}});

    return ($num_nodes_used);
}



####
sub _get_exit_values {
    my $self = shift;
    my $num_cmds = $self->{num_cmds};
    my @retValues;

    #print "Processing $retvals_dir\n";
    for (my $i = 0; $i < $num_cmds; $i++) {

		my $retval_file = $self->_get_ret_filename($i);

        if (-s $retval_file) {
            open (my $fh, $retval_file) or die $!;
            my $retval_string = <$fh>;
            $retval_string =~ s/\s//g;
            $retValues[$i] = $retval_string;
            close $fh;
        } else {
            $retValues[$i] = "FILE_NOT_EXISTS";
        }
    }
    $self->{retvalues} = \@retValues;
}


sub get_failed_cmds {
    my $self = shift;
    my $retvalues_aref = $self->{retvalues};
    my $cmds_list_aref = $self->{cmds_list};

    my @failed_cmds;
    for (my $i = 0; $i <= $#$retvalues_aref; $i++) {
        my $retval = $retvalues_aref->[$i];
        if ($retval) {
            push (@failed_cmds,
                  { cmd => $cmds_list_aref->[$i],
                    ret => $retval,
                } );
        }
    }
    return (@failed_cmds);
}


sub _wait_for_completions {
    my $self = shift;

    print "sub _wait_for_completions()\n" if $self->verbose;

    my $nodes_in_progress_href = $self->{nodes_in_progress};
    my $seen_finished = 0;
    my @done;

    while (! $seen_finished) {

        ## check for finished jobs
        foreach my $monitor_file (keys %$nodes_in_progress_href) {
            if (-e $monitor_file) {
                push (@done, $monitor_file);
                $seen_finished = 1;
            }
        }
        if ($seen_finished) {
            foreach my $monitor_file (@done) {
                my $job_id = $nodes_in_progress_href->{$monitor_file};
                delete $nodes_in_progress_href->{$monitor_file}; #remove from queue
                delete $self->{job_id_to_cmd_indices}->{$job_id};
                delete $self->{job_id_to_submission_time}->{$job_id};
            }
            return (scalar (@done)); #num jobs completed
        }
        else {
            ## wait a while and check again
            sleep(1);
        }
    }
}


sub _write_pid_file {
    my $self = shift;
    my $log_dir = $self->{log_dir};
	  my $hostname = `hostname`;
	  chomp($hostname);
    open (my $fh, ">$log_dir/$hostname.pid") or die $!;
    print $fh $$;
    close $fh;
}


sub _write_result_summary {
    my ($self, $num_successes, $num_failures, $num_unknown) = @_;
    my $status = ($num_failures == 0 && $num_unknown == 0) ? "success" : "failure";

    $self->{status} = $status;
    $self->{num_failures} = $num_failures;
    $self->{num_successes} = $num_successes;
    $self->{num_unknown} = $num_unknown;

    my $log_dir = $self->{log_dir};
    open (my $fh, ">$log_dir/job.finished.$status") or die $!;
    print $fh "num_successes: $num_successes\n"
        . "num_failures: $num_failures\n"
        . "num_unknown: $num_unknown\n";
    close $fh;

}

sub clean_logs {
    my $self = shift;
    my $log_dir = $self->{log_dir};

    my $cmd = "rm -rf $log_dir";
    system $cmd;
    return ($?);
}


sub _write_minimal_environment {
    my ($self, $fh) = @_;

    print $fh <<_EOFENV_;

## add any special environment settings

echo HOST: \$HOSTNAME
echo HOST: \$HOSTNAME >&2

_EOFENV_

;

    return;

}


####
sub _get_ret_filename {
    my $self = shift;
    my ($cmd_index) = @_;

    my $retvals_dir = $self->{retvals_dir};

    my $retval_bin = int ($cmd_index / $self->retval_bin_size);
    my $retval_file = $retvals_dir . "/$retval_bin/entry_$cmd_index.ret";

    return($retval_file);
}


sub run {
    my $self = shift;

    #----------- check if no command to run no need to go further #-----------
    if ($self->num_cmds == 0){
      print "No command to submit to the grid. Quiting now.\n\n";
      return();
    }

    $self->_write_pid_file();

    my $max_job_at_a_time = $self->{max_job_at_a_time};
    my $num_cmds = $self->{num_cmds};
    my $num_cmds_launched = 0;
    my $num_running_job = 0;


    #----------- Submit all jobs -----------
    if($num_cmds > $max_job_at_a_time){
        print STDERR "$max_job_at_a_time max job send to the scheduler at a time. I Will wait that somes jobs are finished before running the nexts.";
    }

    while ($num_cmds_launched < $num_cmds) {
        $num_cmds_launched = $self->_submit_job($num_cmds_launched);
        print STDERR "\r  CMDS: $num_cmds_launched / $num_cmds submitted";
        $num_running_job = $self->_get_num_nodes_used();
        if ($num_running_job >= $max_job_at_a_time) {
            my $num_jobs_finished = $self->_wait_for_completions();
            $num_running_job -= $num_jobs_finished;
        }
    }
    print STDERR "\n* All cmds submitted to grid. Now waiting for them to finish.\n";
    print STDERR "\r  CMDS: $num_cmds_launched / $num_cmds remaining to complete";


    #----------- Wait all job are finished -----------
    while (my $num_jobs_finished = $self->_wait_for_completions()) {
        $num_running_job -= $num_jobs_finished;
        print STDERR "\r  CMDS: $num_running_job / $num_cmds remaining to complete";
        if ($num_running_job == 0){
           last;
        }
    }
    print STDERR "\n* All nodes completed. Now auditing job completion status values\n";


    # ----------- Check the jobs --------------
    $self->_get_exit_values();

    my $retvals_aref = $self->{retvalues};
    my $num_successes = 0;
    my $num_failures = 0;
    my $num_unknown = 0;

    foreach my $retval (@$retvals_aref) {

        if ($retval =~ /\d+/) {
            if ($retval == 0) {
                $num_successes++;
            }
            else {
                $num_failures++;
            }
        }
        else {
          print "retval= $retval\n";
            $num_unknown++;
        }
    }

    $self->_write_result_summary($num_successes, $num_failures, $num_unknown);

    if ($num_successes == $num_cmds) {
        $self->clean_logs();
        print "All $num_cmds completed successfully.\n\n";
    }
    else {
        print "Failures encountered:\n"
            . "num_success: $num_successes\tnum_fail: $num_failures\tnum_unknown: $num_unknown\n";

        my @failed_cmds = $self->get_failed_cmds();
        my $num_failed_cmds = scalar (@failed_cmds);
        my $msg = "$num_failed_cmds of $num_cmds failed = " . ($num_failed_cmds / $num_cmds * 100) . " % failure.\n";

        print STDERR $msg;

        my @failed_cmd_lines;
        open (my $failed_fh, ">failed_cmds.$$") or die $!;
        foreach my $failed_cmd (@failed_cmds) {
            my $cmd = $failed_cmd->{cmd};
            my $ret = $failed_cmd->{ret};
            print $failed_fh "$cmd\nRET($ret)\n\n";
            push (@failed_cmd_lines, $cmd);
        }
        close $failed_fh;
        $self->clean_logs() unless $keep_logs_on_failure;
        print "View file \'failed_cmds.$$\' for list of commands that failed.\n\n";
        return (@failed_cmd_lines); # at least one job failed.
    }

    print "Finished.\n\n";
}

package main;
use strict;
use warnings;
use File::Basename;


#### Test routine, run by executing module directly:
if (basename($0) eq 'Job.pm') {
    my $ret = &Job::run("ls");
    if ($ret) {
        die "Error, test failed.\n";
    }
    else {
        print "Test ran successfully.\n";
    }

    exit($ret);

}

;
1
