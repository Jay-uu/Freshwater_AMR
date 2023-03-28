#!/bin/bash -l
#SBATCH -A naiss2023-5-97
#SBATCH -p core -n 2
#SBATCH -t 8-00:00:00 
#SBATCH -J sourmash
#SBATCH -o /proj/fume/nobackup/private/jay/Freshwater_AMR/scripts/process_sites/logs/sourmash_230327.log 
#SBATCH -e /proj/fume/nobackup/private/jay/Freshwater_AMR/scripts/process_sites/logs/sourmash_230327.err 
#SBATCH --mail-user jay.hakansson.4449@student.uu.se
#SBATCH --mail-type=FAIL,END

READ_DIRS='/proj/fume/raw_data/SITES/SITES_metagenomes_Jan2023'
RESULTS='/proj/fume/private/jay/process_sites/01_sourmash'

module load conda
source conda_init.sh
export CONDA_ENVS_PATH=/proj/fume/nobackup/private/jay/Freshwater_AMR/conda_envs
bash

mamba activate sourmash-4.6.1

cd $TMPDIR
mkdir reads
for dir in $READ_DIRS/Sample_*
do
cp $dir/*.fastq.gz reads/
done

#trying to shorten the names
#testing if suffix is necessary
cd reads
mmv "VK-3496-*" "#1"
mmv "*_L002*" "#1#2"
mmv "*_001.fastq.gz*" "#1#2"
cd ..

mkdir sig_files
sourmash sketch dna reads/*.fastq.gz --outdir sig_files

sourmash compare sig_files/*.sig -o distances.cmp -k 31

sourmash plot distances.cmp --output-dir $RESULTS
cp distances.cmp $RESULTS/
cp -r sig_files $RESULTS/
cp distances.cmp.labels.txt $RESULTS/

