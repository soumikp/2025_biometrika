#!/bin/bash

#SBATCH --mail-type=FAIL
#SBATCH --time=0-02:00
#SBATCH --job-name=collider_comparison  
#SBATCH --mem-per-cpu=2500MB
#SBATCH --array=1-324
#SBATCH --cpus-per-task=1
#SBATCH --output=/ihome/spurkayastha/soumik/2025_biometrika/log/simulation/slurm-%x_%A_%a.out
#SBATCH --account=spurkayastha
   
module load R/4.4.0
Rscript /ihome/spurkayastha/soumik/2025_biometrika/code/simulations.R