#!/usr/bin/env nextflow
/*
========================================================================================
Multiple co-assemblies of metagenomic soil and water samples
========================================================================================
Trimmed reads have been sorted into different directories for different assembly strategies
----------------------------------------------------------------------------------------
*/

process assemble_water {
    publishDir "${params.publish_dir}/03_megahit", mode: 'copy'
    conda "${params.envs_dir}/megahit-1.2.9"
    input:
    path sample
    output:
    tuple val("${sample.baseName}"), path("${sample.baseName}_asm/${sample.baseName}_megahit.fa.gz"), path("*trimmed.fastq.gz")
    """
    cp $sample/*trimmed.fastq.gz .
    sample_ID="${sample.baseName}"
    read1=`ls *R1* | grep -v / | xargs echo | sed 's/ /,/g'`
    read2=`ls *R2* | grep -v / | xargs echo | sed 's/ /,/g'`
    megahit -1 \$read1 -2 \$read2 -o \${sample_ID}_asm --out-prefix \${sample_ID}_megahit --min-contig-len 1000 -t 20
    mv \${sample_ID}_asm/\${sample_ID}_megahit.contigs.fa \${sample_ID}_asm/\${sample_ID}_megahit.fa
    gzip \${sample_ID}_asm/\${sample_ID}_megahit.fa
    """

}

/*
process assemble_soil {
publishDir "${params.publish_dir}/03_megahit", mode: 'copy'
    conda "${params.envs_dir}/megahit-1.2.9"
    input:
    path sample
    output:
    tuple val("${sample.baseName}"), path("${sample.baseName}_megahit.fa.gz"), path("*trimmed.fastq.gz")
    """
    cp $sample/*trimmed.fastq.gz .
    sample_ID="${sample.baseName}"
    read1=`ls *R1* | grep -v / | xargs echo | sed 's/ /,/g'`
    read2=`ls *R2* | grep -v / | xargs echo | sed 's/ /,/g'`
    megahit -1 \$read1 -2 \$read2 -o \${sample_ID}_asm --out-prefix \${sample_ID}_megahit --min-contig-len 1000 -t 20 --presets meta-large
    mv \${sample_ID}_asm/\${sample_ID}_megahit.contigs.fa \${sample_ID}_asm/\${sample_ID}_megahit.fa
    gzip \${sample_ID}_asm/\${sample_ID}_megahit.fa
    """
}
*/

process map_reads_asm {
    publishDir "${params.publish_dir}/04_mapping/$sample/", mode: 'copy', pattern: "*sorted.bam"
    input:
    tuple val(sample), path(asm), path(read_files)
    output:
    tuple val("$sample"),path("$asm"), path("*_map_sorted.bam")
    """
    module load bioinfo-tools bwa-mem2 samtools
    bwa-mem2 index -p $sample $asm

    r1=(*R1*trimmed.fastq.gz)
    r2=(*R2*trimmed.fastq.gz)
    for ((i = 0; i < \${#r1[@]} && i < \${#r2[@]}; i++))
    do
    read_samp=`basename \${r1[i]} _R1_trimmed.fastq.gz`
    bwa-mem2 mem -t 20 $sample "\${r1[i]}" "\${r2[i]}" | samtools view -@ 20 -bo \${read_samp}_map.bam
    samtools sort \${read_samp}_map.bam -o \${read_samp}_map_sorted.bam
    done
    """
}


/*
process binning {
    publishDir "${params.publish_dir}/05_metabat", mode: 'copy'
    input:
    tuple val(sample), path(asm), path(mapped_reads)
    output:
    tuple val("$sample"), path("${sample}_bins")
    """
    module load bioinfo-tools MetaBat/2.12.1
    jgi_summarize_bam_contig_depths --outputDepth ${sample}_depth.txt *_sorted.bam
    mkdir ${sample}_bins
    runMetaBat.sh -i ${asm} -a ${sample}_depth.txt -o ${sample}_bins/${sample}_ --unbinned
    """
}

process rename_contigs {
    publishDir "${params.publish_dir}/05_metabat", mode: 'copy'
    input:
    tuple val(sample), path(bins_dir)
    output:
    tuple val("$sample"), path()
    """
    //#!/usr/bin/env python
    from Bio import SeqIO
    import math
    import os
    for bin in os.listdir($bins_dir):
    
    """
}
*/

workflow {
    params.publish_dir="/proj/fume/private/jay/process_sites"
    params.envs_dir="/proj/fume/nobackup/private/jay/Freshwater_AMR/conda_envs"
    //groups are directories with multiple trimmed reads, while singles are trimmed PE reads for individual samples
    params.water_groups="/crex/proj/fume/private/jay/process_sites/02_fastp/group_sets/water"
    //params.soil_groups="/crex/proj/fume/private/jay/process_sites/02_fastp/group_sets/soil"
    params.water_singles="/crex/proj/fume/private/jay/process_sites/02_fastp/water"
    //params.soil_singles="/crex/proj/fume/private/jay/process_sites/02_fastp/soil"
    //for testing
    //water = Channel.fromPath( "${params.water_groups}/Abisko", type:'dir' )
    water = Channel.fromPath( "${params.water_singles}/Sample-A1/", type:'dir'  )
    assemble_water(water) | map_reads_asm  
}

