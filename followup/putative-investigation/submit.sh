# this script is to generate and submit job array in slurm and dSQ system
# we recommend you cd into a directory for storing job file and dSQ batch script
# then copy paste the content of this script 
# do not run this script directly

module load dSQ


partition="pi_hall"

nthreads=12
memory=60G


experiment_name="3-putative-investigation"
# this is the filename of parameter file

job_file_name="3-putative-investigation_batch.txt"
# this is the jobe file name you want to generate

grand_dir="/home/xw445/scripts/LOCATER_real_data_vignette/"

log_dir="/home/xw445/palmer_scratch/logs/WashU_CCDG/to_share/followup/${experiment_name}/"
mkdir -p $log_dir
# please specify your own log directory


######## generating job file #############
file_path="/home/xw445/gibbs/WashU_CCDG/screening/to_share/locater_screening_local_data/wg_unrelated_METSIMonlyData_101phenos_Bjitter/sig_locater_candidate_loci_all_clean.txt"
# file that contain chr, start and end positions of putative loci

tail -n +2 "$file_path" | while read -r line; do
  # skip header, Extract columns into variables
  chr=$(echo "$line" | awk '{print $1}')
  start=$(echo "$line" | awk '{print $2}')
  end=$(echo "$line" | awk '{print $3}')

  echo "export APPTAINER_BIND="/home/xw445:/home/xw445"; apptainer exec /home/xw445/scripts/container-setup/mini-shark-nov28 ${grand_dir}followup/putative-investigation/run_jobs followup ${experiment_name} ${nthreads} ${chr} ${start} ${end}" >> ${job_file_name}

done

# /home/xw445/scripts/container-setup/mini-shark-nov28: the downloaded docker container as Apptainer image
# (https://docs.ycrc.yale.edu/clusters-at-yale/guides/containers/) 
# ${grand_dir}whole_genome_screening/putative-investigation/run_jobs: executable file that runs a R script with some environment variables
# ${experiment_name}: filename (no extensions) of R script file

######## finish generating job file #############


######## submit job array  #############

dsq --job-file ${job_file_name} --cpus-per-task ${nthreads} --mem ${memory} -t 96:00:00 \
--mail-type ALL --status-dir ${log_dir} --max-jobs 45 \
--partition ${partition} --nodes 1 --output ${log_dir}%A_%a.out
# specify resource usage as needed

######## finish submitting job array  #############


