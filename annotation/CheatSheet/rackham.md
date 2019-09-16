# Rackham

Rackham is a high performance computer cluster at UPPMAX (Uppsala Multidisciplinary Center for Advanced Computational Science).
UPPMAX is Uppsala University's resource of high-performance computers, large-scale storage and know-how of high-performance computing (HPC).

  * [UPPMAX github repo](https://github.com/UPPMAX?utf8=✓&q=&type=&language=)
  * [UPPMAX official Cheat Sheet](https://www.uppmax.uu.se/support/getting-started/uppmax-cheat-sheet/)
  * [Rackham User Guide](https://uppmax.uu.se/support/user-guides/rackham-user-guide/)
  

### UPPMAX modules

| Command | Comment
| --- | --- 
| module avail   | List available modules
| module load ***modulename***   | Load the module ***modulename***
| module unload ***modulename***   | Load the module ***modulename***
| module list  | List all modules loaded
| module spider ***sam***  | Search all module containing ***sam***

### Running jobs with the Slurm ressource manager

| Command | Comment
| --- | --- 
| interactive -A project   | Start interactive job
| sbatch -A projectID -t d-hh:mm:ss -n cores -p partition my_jobscript_file  | Start batch job
| sbatch -A projectID -t 7-00:00:00 -n 16 -p node my_jobscript_file  | Running for 7 days on 16 cores one node partition
| scancel ***jobId***   | Cancel a single job
| scancel -i -u user | Interactively cancel all jobs for user

### Showing user, job and project info

Those commans are stored in */sw/uppmax/bin/*

| Command | Comment | Full documentation
| --- | --- | ---
| jobinfo   | Show all running and waiting jobs in the queue | [here](https://www.uppmax.uu.se/support/user-guides/jobstats-user-guide/)
| jobinfo -u ***user***   | Show jobs for specific user | [here](https://www.uppmax.uu.se/support/user-guides/jobstats-user-guide/)
| jobinfo -p -A ***projectID***   | generate plots for resource usage for jobs | [here](https://www.uppmax.uu.se/support/user-guides/jobstats-user-guide/)
| uquota | Show current user's disk usage |
| projinfo    | Show used core hours for current user's projects |
| egrep '^b2011999' /etc/slurm/grantfile | View details of a specific project |
| finishedjobinfo    | Telling you about finished jobs on Rackham |
| projmembers    | Telling you about project membership |
| projsummary  ***projectID***  | Summarizes some useful information about projects |
| jobstats  ***jobID***  | Discover jobstats for the specified job |
| jobstats  -A ***projectID***  | Discover jobstats for the specified project |
| squeue -u ***user*** | Information about ***user***'s jobs |
