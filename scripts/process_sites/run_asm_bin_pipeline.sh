#!/bin/bash -l
#SBATCH -A naiss2023-5-97
#SBATCH -p core -n 2
#SBATCH -t 10-00:00:00 
#SBATCH -J nf_asm_bin
#SBATCH -o /proj/fume/nobackup/private/jay/Freshwater_AMR/scripts/process_sites/logs/nextflow_asm20230510.log
#SBATCH -e /proj/fume/nobackup/private/jay/Freshwater_AMR/scripts/process_sites/logs/nextflow_asm20230510.err
#SBATCH --mail-user jay.hakansson.4449@student.uu.se
#SBATCH --mail-type=FAIL,END

cd /proj/fume/nobackup/private/jay/Freshwater_AMR/scripts/process_sites/
module load conda
source conda_init.sh
export CONDA_ENVS_PATH=/proj/fume/nobackup/private/jay/Freshwater_AMR/conda_envs
bash

mamba activate nextflow-22.10.6
nextflow run asm_bin_pipeline.nf -c asm_bin_pipeline.config -resume
