#!/bin/bash -l
#SBATCH -A naiss2023-5-97
#SBATCH -p core -n 2
#SBATCH -t 4-00:00:00 
#SBATCH -J extract_seqs
#SBATCH -o /proj/fume/nobackup/private/jay/Freshwater_AMR/scripts/amr_analysis/logs/extract_seqs20230421.log 
#SBATCH -e /proj/fume/nobackup/private/jay/Freshwater_AMR/scripts/amr_analysis/logs/extract_seqs20230421.err 
#SBATCH --mail-user jay.hakansson.4449@student.uu.se
#SBATCH --mail-type=FAIL,END

module load bioinfo-tools SeqKit
#v  0.15.0

#cd $TMPDIR
cd /home/jay/master_thesis/results/amr_analysis/
#mkdir tmp_buck_bins
#cd tmp_buck_bins

echo '==================== Start ================'
date

#assign variables
BINS='/home/moritz/proj_folder/uppstore2018116/webexport/stratfreshdb/bins'
RESULT='/home/jay/master_thesis/results/amr_analysis/buck_seqs/'

#file with seq_names like this: ErkenSummer_megahit_metabat_bin-0252_01352
SEQ_IDS='/home/jay/master_thesis/results/amr_analysis/seq_ids_amr_hits_buck.txt'

#file with bin names like this: ErkenSummer_megahit_metabat_bin-3118
ABR_BIN_NAMES='/home/jay/master_thesis/results/amr_analysis/bin_ids.txt'

#mkdir abr_seqs
#for bin in $BINS/*/*
#do
#echo "Extracting bin: $bin"
#tar -xf ${bin} -C . --wildcards '*.faa.gz' --strip-components=3
#done

#cat $ABR_BIN_NAMES | while read file_name
for file_name in `cat $ABR_BIN_NAMES`
do
   echo "Finding sequences for $file_name"
   seqkit grep -f $SEQ_IDS tmp_buck_bins/${file_name}.faa.gz -o buck_seqs/${file_name}.faa.gz
done

#copy res
#cp -r abr_seqs/ $RESULT/

echo '====================== Stop ======================'
date

