# this script is to generate and submit job array in slurm and dSQ system
# we recommend you cd into a directory for storing job file and dSQ batch script
# then copy paste the content of this script 
# do not run this script directly

module load dSQ


partition="pi_hall"

nthreads=4
memory=20G


experiment_name="2-pval-rank-matched"
# this is the filename of parameter file

job_file_name="2-pval-rank-matched_batch.txt"
# this is the job file name you want to generate

grand_dir="/home/xw445/scripts/LOCATER_real_data_vignette/"

log_dir="/home/xw445/palmer_scratch/logs/WashU_CCDG/to_share/rank_matching_selection/${experiment_name}/"
mkdir -p $log_dir
# please specify your own log directory


######## generating job file #############

for i in {1..4000} # our genome was seperated to 4000 small regions
do
	echo "export APPTAINER_BIND="/home/xw445:/home/xw445"; apptainer exec /home/xw445/scripts/container-setup/mini-shark-nov28 ${grand_dir}rank_matching_selection/run_jobs rank_matching_selection ${experiment_name} ${nthreads} ${i}" >> ${job_file_name}

done

# /home/xw445/scripts/container-setup/mini-shark-nov28: the downloaded docker container as Apptainer image
# (https://docs.ycrc.yale.edu/clusters-at-yale/guides/containers/) 
# ${grand_dir}rank_matching_selection/run_jobs: executable file that runs a R script with some environment variables
# ${experiment_name}: filename (no extensions) of R script file

######## finish generating job file #############


######## submit job array  #############

dsq --job-file ${job_file_name} --cpus-per-task ${nthreads} --mem ${memory} -t 24:00:00 \
--mail-type ALL --status-dir ${log_dir} --max-jobs 80 \
--partition ${partition} --nodes 1 --output ${log_dir}%A_%a.out
# specify resource usage as needed

######## finish submitting job array  #############


