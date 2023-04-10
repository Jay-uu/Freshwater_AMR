#!/bin/bash -l
#SBATCH -A naiss2023-5-97
#SBATCH -p core -n 2
#SBATCH -t 5-00:00:00 
#SBATCH -J nextflow_trim
#SBATCH -o /proj/fume/nobackup/private/jay/Freshwater_AMR/scripts/process_sites/logs/nextflow_trim20230410.log #write slurm log name here
#SBATCH -e /proj/fume/nobackup/private/jay/Freshwater_AMR/scripts/amr_finding/logs/nextflow_trim20230410.err #write same slurm log name here
#SBATCH --mail-user jay.hakansson.4449@student.uu.se
#SBATCH --mail-type=FAIL,END

cd /proj/fume/nobackup/private/jay/Freshwater_AMR/scripts/process_sites/
module load conda
source conda_init.sh
export CONDA_ENVS_PATH=/proj/fume/nobackup/private/jay/Freshwater_AMR/conda_envs
bash

mamba activate nextflow-22.10.6
nextflow run trim_pipeline.nf -c trim_pipeline.config -resume
