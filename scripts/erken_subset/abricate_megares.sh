#!/bin/bash -l
echo '==================== Start ================'
date

#assign variables 
WORK_DIR=/home/jay/work_dir
BINS=/home/moritz/data/data_submit/bins/ErkenSummer
RESULT=/home/jay/data/erken_results/02_abricate

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

#activate environment
conda activate abricate_env

echo '======================= Run Abricate ======================='
#run abricate
abricate --db megares --minid 50 bins/*.fna.gz > res/megares50_results.tab

echo '======================= Done running Abricate ==================='
#move results to slow drive
cp -r res/* $RESULT/
rm -r res
#if copied data, delete
rm -r bins

#turn off environment
conda deactivate

echo '====================== Stop ======================'
date
