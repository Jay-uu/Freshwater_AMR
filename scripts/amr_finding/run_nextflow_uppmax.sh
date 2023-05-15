#!/bin/bash -l
#SBATCH -A naiss2023-5-97
#SBATCH -p core -n 2
#SBATCH -t 5-00:00:00 
#SBATCH -J nextflow_amr
#SBATCH -o /proj/fume/nobackup/private/jay/Freshwater_AMR/scripts/amr_finding/logs/nf_amr_pipeline_230512.log
#SBATCH -e /proj/fume/nobackup/private/jay/Freshwater_AMR/scripts/amr_finding/logs/nf_amr_pipeline.err
#SBATCH --mail-user jay.hakansson.4449@student.uu.se
#SBATCH --mail-type=FAIL,END

cd /proj/fume/nobackup/private/jay/Freshwater_AMR/scripts/amr_finding
module load conda
source conda_init.sh
export CONDA_ENVS_PATH=/proj/fume/nobackup/private/jay/Freshwater_AMR/conda_envs
bash

mamba activate nextflow-22.10.6
nextflow run amr_finding_pipeline.nf -c amr_finding_uppmax.config -resume
