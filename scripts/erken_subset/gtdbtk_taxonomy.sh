#!/bin/bash -l

echo '==================== Start ================'
date

#assign variables
WORK_DIR=/home/jay/work_dir
BINS=/home/moritz/data/data_submit/bins/ErkenSummer
RESULT=/home/jay/data/erken_results/04_GTDB-Tk
DB=$1
#go to work directory
cd $WORK_DIR
mkdir res
mkdir bins

echo '======================= Exctracting bins ===================='
#extract all bins fasta files: .fna.gz
for dir in $BINS/*
do
tar -xf $dir -C bins --wildcards '*.fna.gz' --strip-components=3
done

date
#activate environment
conda activate gtdbtk_env

echo '======================= Run GTDB-Tk ======================='
#add command


echo '======================= Done running GTDB-Tk ==================='
#move results to slow drive
cp res/* $RESULT/
#rm -r res #do this manually
#delete copied data
rm -r bins

#turn off environment
conda deactivate

echo '====================== Stop ======================'
date
