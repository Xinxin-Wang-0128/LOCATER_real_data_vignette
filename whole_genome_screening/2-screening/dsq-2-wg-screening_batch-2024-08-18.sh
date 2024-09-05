#!/bin/bash
#SBATCH --output /home/xw445/palmer_scratch/logs/WashU_CCDG/to_share/whole_genome_screening/2-wg-screening/%A_%a.out
#SBATCH --array 0-3%45
#SBATCH --job-name dsq-2-wg-screening_batch
#SBATCH --cpus-per-task 12 --mem 60G -t 24:00:00 --mail-type ALL --partition pi_hall --nodes 1

# DO NOT EDIT LINE BELOW
/vast/palmer/apps/avx2/software/dSQ/1.05/dSQBatch.py --job-file /vast/palmer/home.mccleary/xw445/scripts/shark/WashU_CCDG/to_share/whole_genome_screening/2-screening/2-wg-screening_batch.txt --status-dir /home/xw445/palmer_scratch/logs/WashU_CCDG/to_share/whole_genome_screening/2-wg-screening

