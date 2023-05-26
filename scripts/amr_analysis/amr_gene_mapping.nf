#!/usr/bin/env nextflow

#script by Moritz Buck
#https://github.com/moritzbuck


nextflow.enable.dsl=2

params.database = "/home/moritz/projects/0073_jays_maps/data/clustered_amrs"
params.data_folder = "/home/moritz/proj_folder/uppstore2017149/webexport/stratfreshdb/reads/mg/"
params.running_folder = "/home/moritz/projects/0073_jays_maps"

workflow {
    Channel.fromFilePairs("${params.data_folder}/*_{fwd,rev}.fastq.gz", checkIfExists:true) | GENEMAPPER
  }

process GENEMAPPER {
    debug true

    input:
     tuple val(s), file(reads)


    output:
    file("result/${s}.mmseqs2_out")

    script:
    """
    cd ${params.running_folder}
    echo creating $s db
    mkdir -p  $SNIC_TMP
    mmseqs createdb ${params.data_folder}/${s}_fwd.fastq.gz ${params.data_folder}/${s}_rev.fastq.gz $SNIC_TMP/${s} 2>&1 >> logs/${s}
    echo mapping ${s} 
    mmseqs map --threads 20  $SNIC_TMP/${s} ${params.database} $SNIC_TMP/${s}.map $SNIC_TMP  2>&1 >> logs/${s}
    echo converting output to blastout
    mmseqs convertalis  $SNIC_TMP/$s data/clustered_amrs $SNIC_TMP/${s}.map result/${s}.mmseqs2_out  --threads 20  2>&1  >> logs/${s}
    rm -r $SNIC_TMP/${s}*

    """
}
