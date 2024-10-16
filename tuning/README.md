This folder contains example scripts for tuning. 
To start tuning, you will need to start from submit.sh.

submit.sh: This Bash script is designed to automate the process of generating and submitting a job array using a Slurm workload manager with dSQ, a tool for running job arrays on Slurm clusters.
(more information: https://docs.ycrc.yale.edu/clusters-at-yale/job-scheduling/dsq/)

In submit.sh, we generate a job file, then submit jobs in this job file in a job array. 
each line of the job file refers to a job, which in our case will use to our tuning parameter R script and a executable file.
(submit.sh -> run_tuning -> R script in tuning_pars)

run_tuning: the executable file that sets environment variables that limit the number of threads used by R and OpenBLAS,
Constructs the path to the R script based on input arguments and then executes the R script using the Rscript command with additional arguments.
(R script is the tuning parameter residing in tuning_pars)

tuning_data_collect.R is used for collecting the tuning data.
tuning_viz.R is used for visualizing the tuning data. We recommend using it in Rstudio since we would like to generate interactive plots saved in HTML format. It is easier to make standalone HTML files in Rstudio.
