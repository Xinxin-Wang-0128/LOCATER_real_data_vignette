This folder contains example scripts for whole genome screening and putative loci investigation.


2. screening
To start, you will need to start from submit.sh.

submit.sh: This Bash script is designed to automate the process of generating and submitting a job array using a Slurm workload manager with dSQ, a tool for running job arrays on Slurm clusters.
(more information: https://docs.ycrc.yale.edu/clusters-at-yale/job-scheduling/dsq/)

In submit.sh, we generate a job file, then submit jobs in this job file in a job array. 
each line of the job file refers to a job, which in our case will use to our screening R script and a executable file.
(submit.sh -> run_jobs in 2-screening folder -> R script for screening in 2-screening folder)

run_jobs: the executable file that sets environment variables that limit the number of threads used by R and OpenBLAS,
Constructs the path to the R script based on input arguments and then executes the R script using the Rscript command with additional arguments.
(R script is the screening script residing in 2-screening folder)

3. putative investigation
To start, you will need to start from submit.sh.

In submit.sh, we generate a job file, then submit jobs in this job file in a job array. 
each line of the job file refers to a job, which in our case will use to our screening R script and a executable file.
(submit.sh -> run_jobs in 3-putative-investigation folder -> R script for screening in 3-putative-investigation folder)



