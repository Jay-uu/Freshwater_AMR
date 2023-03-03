#!/bin/bash -l
#write the database name you want to use on the command line when you start the script
echo '==================== Start ================'
date

#assign variables
WORK_DIR='/home/jay/work_dir'
BINS='/home/moritz/data/data_submit/bins/ErkenSummer'
RESULT='/home/jay/data/erken_results/04_comp_aro_hits'
SEQ_IDS='/home/jay/data/erken_results/04_comp_aro_hits/abr_seq_ids.txt'
ABR_BIN_NAMES='/home/jay/data/erken_results/04_comp_aro_hits/abr_bin_names.txt'
#go to work directory
cd $WORK_DIR
mkdir abr_seqs

mamba activate seqkit

cat $ABR_BIN_NAMES | while read file_name 
do
   BIN_ID=`basename $file_name .fna.gz`
   echo "Extracting bin: $BIN_ID" 
   tar -xf ${BINS}/${BIN_ID}.tar.gz -C . --wildcards '*.faa.gz' --strip-components=3
   echo "Finding sequences"
#this part not working. Result files are empty.
   seqkit grep -f $SEQ_IDS ${BIN_ID}.faa.gz -o abr_seqs/${BIN_ID}.fa.gz
done

#move results to slow drive
cp -r abr_seqs $RESULT/
#clean work dir
rm -r abr_seqs
rm -r *.faa.gz

#turn off environment
mamba deactivate

echo '====================== Stop ======================'
date
