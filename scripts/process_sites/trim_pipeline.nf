#!/usr/bin/env nextflow
/*
========================================================================================
Trimming of raw reads
========================================================================================
----------------------------------------------------------------------------------------
*/


process fastp {
    publishDir "${params.publish_dir}/02_fastp", mode: 'copy'
    conda "${params.envs_dir}/fastp-0.23.2"
    input:
    path dir
    output:
    tuple val("${dir.baseName}"), path("${dir.baseName}_R1_trimmed.fastq.gz"), path("${dir.baseName}_R2_trimmed.fastq.gz")
    """
    fastp --in1 $dir/*R1*.fastq.gz --in2 $dir/*R2*.fastq.gz --out1 ${dir.baseName}_R1_trimmed.fastq.gz --out2 ${dir.baseName}_R2_trimmed.fastq.gz --detect_adapter_for_pe --trim_poly_g --low_complexity_filter
    """

}


workflow {
    params.publish_dir="/proj/fume/private/jay/process_sites"
    params.envs_dir="/proj/fume/nobackup/private/jay/Freshwater_AMR/conda_envs"
    //params.read_dirs="/proj/fume/raw_data/SITES/SITES_metagenomes_Jan2023"

    raw_reads = Channel.fromPath( "/proj/fume/raw_data/SITES/SITES_metagenomes_Jan2023/*" , type: 'dir')

    fastp(raw_reads)
}
