#!/bin/bash -l
date
#assign variables
WORK_DIR=/home/jay/work_dir
BINS=/home/moritz/data/data_submit/bins/ErkenSummer
RESULT=/home/jay/data/erken_results/03_rgi

#go to work directory
cd $WORK_DIR
mkdir res
mkdir bins

echo '=========== Extracting files ========='
#extract all bins fasta files: .fna.gz
for dir in $BINS/*
do
tar -xf $dir -C bins --wildcards '*.fna.gz' --strip-components=3
done

#activate environment
conda activate rgi_env

echo '========== Run RGI =========='
for file in bins/*
do
   BIN_ID= basename $file .fna.gz
   echo 'running bin:'
   echo "$BIN_ID"
   rgi main --input_sequence $file --output_file res/$BIN_ID.tab -n 23 --include_loose -a DIAMOND --low_quality
done

echo '========== RGI done =========='
#move results to slow drive
cp -r res $RESULT/
#rm -r res
#delete copied data
rm -r bins

#turn off environment
conda deactivate

date
