#!/bin/bash -l

#assign variables
WORK_DIR=/home/jay/work_dir
BINS=/home/moritz/data/data_submit/bins/ErkenSummer
RESULT=/home/jay/data/erken_results/01_checkm

#go to work directory
cd $WORK_DIR
mkdir res
mkdir bins

#extract all bins fasta files: .fna.gz
for dir in $BINS/*
do
tar -xf $dir -C bins --wildcards '*.fna.gz' --strip-components=3
done
#unzip
gunzip bins/*.fna.gz

#load environment
conda activate checkm_env 

#run checkM
checkm lineage_wf -t 20  bins res > summary.txt

#move results to slow drive 
cp -r res $RESULT/
#rm -r res #or remove manually?
#if copied data, delete
rm -r bins

#turn off environment
conda deactivate
