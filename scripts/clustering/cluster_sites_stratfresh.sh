#!/bin/bash -l
#SBATCH -A naiss2023-5-97
#SBATCH -p core -n 20
#SBATCH -t 7-00:00:00 
#SBATCH -J clustering
#SBATCH -o /proj/fume/nobackup/private/jay/Freshwater_AMR/scripts/clustering/logs/clustering_sites_stratfresh_230525.log
#SBATCH -e /proj/fume/nobackup/private/jay/Freshwater_AMR/scripts/clustering/logs/clustering_sites_stratfresh_230525.err
#SBATCH --mail-user jay.hakansson.4449@student.uu.se
#SBATCH --mail-type=FAIL,END

module load conda
source conda_init.sh
export CONDA_ENVS_PATH=/proj/fume/nobackup/private/jay/Freshwater_AMR/conda_envs
bash

mamba activate mOTUlizer-0.3.2

STRATFRESH_BINS='/crex/proj/uppstore2017149/webexport/stratfreshdb/bins'
STRATFRESH_CHECKM='/crex/proj/fume/private/jay/amr_finding/01_checkm'

SITES_BINS='/crex/proj/fume/private/jay/process_sites/07_bins'
SITES_CHECKM='/crex/proj/fume/private/jay/process_sites/08_amr_finding/01_checkm'

RESULT='/crex/proj/fume/private/jay/clustering'


cd $RESULT

#extract all bins
#mkdir tmp_bins

#for bin in $SITES_BINS/*/*
#do
#echo "Extracting bin: $bin"
#tar -xf ${bin} -C tmp_bins --wildcards '*.fna.gz' --strip-components=2
#done
#echo "=====Done extracting SITES====="

#for bin in $STRATFRESH_BINS/*/*
#do
#echo "Extracting bin: $bin"
#tar -xf ${bin} -C tmp_bins --wildcards '*.fna.gz' --strip-components=3
#done

#echo "=====Done extracting stratfresh====="
#echo "unzipping bins"
#gunzip -r tmp_bins
#printf '%s\n' tmp_bins/*.fna > bin_names.txt

#echo "concatenating checkm-files"
#cat $STRATFRESH_CHECKM/*.txt >> checkm.txt
#cat $SITES_CHECKM/*.txt >> checkm.txt

echo "===== Running motulizer ====="

#run motulizer
mOTUlize.py --fnas bin_names.txt --output mOTUs_sites_stratfresh.tsv --checkm checkm.txt --keep-simi-file sites_stratfresh_simi.txt --threads 20 --txt

#rm -r tmp_bins
#rm checkm.txt

echo "=====Done====="

