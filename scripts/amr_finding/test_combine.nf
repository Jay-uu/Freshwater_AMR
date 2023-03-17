#!/usr/bin/env nextflow
/*
========================================================================================
AMR-finding in multiple metagenomic assembled lake samples
========================================================================================
With quality eval and taxonomic identification
----------------------------------------------------------------------------------------
*/

process checkm_process{
    publishDir "${params.publish_dir}/test1", mode: 'copy'
    input:
    path bin
    output:
    //tuple val("bin_ID"), path ("bin_list.txt")
    //tuple stdout, path ("bin_list.txt")
    tuple val("${bin.baseName}"), path ("bin_list.txt")
    """
    mkdir test
    echo $bin > bin_list.txt
    for b in ${bin}/*.tar.gz
    do
        echo \$b >> bin_list.txt
    done
    #echo "testing"
    bin_ID="${bin.baseName}"
    #this method doesn't work if *anything* other than the bin_id is printed to stdout...
    #printf "\$bin_ID"
    """
}

process abricate_process{
    publishDir "${params.publish_dir}/test2", mode: 'copy'
    input:
    path bin
    output:
    tuple stdout, path ("only_bin.txt")
    //tuple val("bin_ID"), path ("only_bin.txt")
    """
    echo $bin > only_bin.txt
    bin_ID="${bin.baseName}"
    printf "\$bin_ID"
    """
}

workflow {

    Channel
        .fromPath("/home/moritz/data/data_submit/bins/ErkenS*", type:'dir')
        .multiMap { bin_dirs -> bin_dirs_checkm: bin_dirs_abricate: bin_dirs }
        .set { dirs }
    //results directory
    params.publish_dir='/home/jay/data/amr_finding'
    
    checkm_process(dirs.bin_dirs_checkm)
    abricate_process(dirs.bin_dirs_abricate)
    
    checkm_process.out \
        | combine( abricate_process.out, by: 0)
        | view()
}
