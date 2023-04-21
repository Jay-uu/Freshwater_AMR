#!/bin/bash -l
#SBATCH -A naiss2023-5-97
#SBATCH -p core -n 20
#SBATCH -t 10-00:00:00
#SBATCH -J abisko_asm
#SBATCH -o /proj/fume/nobackup/private/jay/Freshwater_AMR/scripts/process_sites/logs/abisko_asm_230419.log
#SBATCH -e /proj/fume/nobackup/private/jay/Freshwater_AMR/scripts/process_sites/logs/abisko_asm_230419.err
#SBATCH --mail-user jay.hakansson.4449@student.uu.se
#SBATCH --mail-type=FAIL,END

READS="/proj/fume/private/jay/process_sites/02_fastp"
OUTDIR="/proj/fume/private/jay/process_sites/03_assemblies"

module load bioinfo-tools megahit

#Samples that are part of set: A1, A3, A4, A5, A6, A7. All are lake.
#Meta-file with set-info: /proj/fume/metadata/sites.csv
#MEGAHIT v1.2.9
cd $TMPDIR

for s in 'A1' 'A3' 'A4' 'A5' 'A6' 'A7'
do
  cp $READS/*${s}*_trimmed.fastq.gz $TMPDIR
done

megahit -1 Sample_VK-3496-A1_R1_trimmed.fastq.gz -2 Sample_VK-3496-A1_R2_trimmed.fastq.gz \
-1 Sample_VK-3496-A3_R1_trimmed.fastq.gz -2 Sample_VK-3496-A3_R2_trimmed.fastq.gz \
-1 Sample_VK-3496-A4_R1_trimmed.fastq.gz -2 Sample_VK-3496-A4_R2_trimmed.fastq.gz \
-1 Sample_VK-3496-A5_R1_trimmed.fastq.gz -2 Sample_VK-3496-A5_R2_trimmed.fastq.gz \
-1 Sample_VK-3496-A6_R1_trimmed.fastq.gz -2 Sample_VK-3496-A6_R2_trimmed.fastq.gz \
-1 Sample_VK-3496-A7_R1_trimmed.fastq.gz -2 Sample_VK-3496-A7_R2_trimmed.fastq.gz \
-o abisko_asm --out-prefix abisko_megahit -t 20

rm *_trimmed.fastq.gz

cp -r * $OUTDIR
