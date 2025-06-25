#!/bin/bash

#SBATCH --job-name=vus_rand
#SBATCH --partition=cpuq
#SBATCH --time=300:00:00
#SBATCH --nodes=1
#SBATCH --output=vus_rand_%a.out
#SBATCH --error=vus_rand_%a.err
#SBATCH --mem-per-cpu=3000MB
#SBATCH --cpus-per-task=30
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aurora.savino@fht.org

module load R/4.1.0-rstudio-dependencies

cd /home/aurora.savino/VUS/VUS/

# run
Rscript ./pipelines/_VUS2024build_00_main_rand.R ${SLURM_ARRAY_TASK_ID}  
