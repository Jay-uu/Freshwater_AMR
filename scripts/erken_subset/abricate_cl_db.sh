#!/bin/bash -l
#write the database name you want to use on the command line when you start the script
echo '==================== Start ================'
date

#assign variables 
WORK_DIR=/home/jay/work_dir
BINS=/home/moritz/data/data_submit/bins/ErkenSummer
RESULT=/home/jay/data/erken_results/02_abricate
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
conda activate abricate_env

echo '======================= Run Abricate ======================='
#run abricate
abricate --db $DB bins/*.fna.gz > res/${DB}_results.tab

echo '======================= Done running Abricate ==================='
#move results to slow drive
cp res/* $RESULT/
rm -r res
#delete copied data
rm -r bins

#turn off environment
conda deactivate

echo '====================== Stop ======================'
date
