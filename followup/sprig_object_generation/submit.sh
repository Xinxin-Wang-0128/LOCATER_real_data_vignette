#rm wg_4587segs_new_pass4.txt
module load dSQ

#nthreads=20
#memory=119G
#partition="general"

#nthreads=18
#memory=90G
partition="pi_hall"

nthreads=12
memory=60G


experiment_name="sprig_object_generation"

job_file_name="sprig_object_generation.txt"

grand_dir="/home/xw445/scripts/LOCATER_real_data_vignette/"

log_dir="/home/xw445/palmer_scratch/logs/WashU_CCDG/screening/sprig_investigation/${experiment_name}/"
mkdir -p $log_dir

file_path="/home/xw445/gibbs/WashU_CCDG/screening/data/whole_genome_yale/locater_screening_local_data/wg_unrelated_METSIMonlyData_101phenos_Bjitter_local_10cM_clean/rerun_loci.txt"
# Skip the first line (header)
tail -n +2 "$file_path" | while read -r line; do
  # Extract columns into variables
  chr=$(echo "$line" | awk '{print $1}')
  start=$(echo "$line" | awk '{print $2}')
  end=$(echo "$line" | awk '{print $3}')
  pos=$(echo "$line" | awk '{print $4}')

 echo "export APPTAINER_BIND="/home/xw445:/home/xw445"; apptainer exec /home/xw445/scripts/container-setup/mini-shark-nov28 ${grand_dir}followup/sprig_object_generation/run_sprig_object_generation followup ${experiment_name} ${nthreads} ${chr} ${start} ${end} ${pos}" \
  >> ${job_file_name}

done


### only edit the following part ####

dsq --job-file ${job_file_name} --cpus-per-task ${nthreads} --mem ${memory} -t 1-00:00:00 \
--mail-type ALL --status-dir ${log_dir} --max-jobs 9 \
--partition ${partition} --nodes 1 --output ${log_dir}%A_%a.out
 
#--max-jobs 34 in pi_hall
# this file is to generate the txt job file, which is a input for gen-dsq-joblist.sh
# need to copy and paste and just run, don't bash this script