#!/bin/bash -l

echo '==================== Start ================'
date

#assign variables
WORK_DIR=/home/jay/work_dir
BINS=/home/moritz/data/data_submit/bins/ErkenSummer
RESULT=/home/jay/data/erken_results/05_taxonomy
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
mamba activate gtdbtk-2.1.1

echo '======================= Run GTDB-Tk ======================='
#add command
gtdbtk classify_wf --genome_dir bins --out_dir res --cpus 20 -x gz

echo '======================= Done running GTDB-Tk ==================='
#move results to slow drive
cp -r res/* $RESULT/
#rm -r res #do this manually
#delete copied data
rm -r bins

#turn off environment
mamba deactivate

echo '====================== Stop ======================'
date
