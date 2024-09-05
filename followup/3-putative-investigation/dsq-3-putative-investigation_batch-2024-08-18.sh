#!/bin/bash
#SBATCH --output /home/xw445/palmer_scratch/logs/WashU_CCDG/to_share/whole_genome_screening/3-putative-investigation/%A_%a.out
#SBATCH --array 0-11%45
#SBATCH --job-name dsq-3-putative-investigation_batch
#SBATCH --cpus-per-task 12 --mem 60G -t 96:00:00 --mail-type ALL --partition pi_hall --nodes 1

# DO NOT EDIT LINE BELOW
/vast/palmer/apps/avx2/software/dSQ/1.05/dSQBatch.py --job-file /vast/palmer/home.mccleary/xw445/scripts/shark/WashU_CCDG/to_share/whole_genome_screening/3-putative-investigation/3-putative-investigation_batch.txt --status-dir /home/xw445/palmer_scratch/logs/WashU_CCDG/to_share/whole_genome_screening/3-putative-investigation

