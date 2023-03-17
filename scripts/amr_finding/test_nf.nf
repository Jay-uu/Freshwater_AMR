#!/usr/bin/env nextflow
/*
========================================================================================
AMR-finding in multiple metagenomic assembled lake samples
========================================================================================
With quality eval and taxonomic identification
----------------------------------------------------------------------------------------
*/


// Get archives with bins
//path to all samples with individual bins in subdirectories per sample
//bin_dirs=Channel.fromPath("/home/moritz/data/data_submit/bins/*", type:'dir')

//start with one bin. I'm thinking ErkenSummer because then I can compare and see if the results are the same!
bin_dirs=Channel.fromPath("/home/moritz/data/data_submit/bins/Erken*", type:'dir')
//split into multiple channels for each tool
bin_dirs.into { bin_dirs_checkm; bin_dirs_abricate }
//results directory
params.publish_dir='/home/jay/data/amr_finding'

process checkm_process{
    publishDir "${params.publish_dir}/test1", mode: 'copy'
    input:
    path bin from bin_dirs_checkm
    output:
    tuple val("bin_ID"), path ("bin_list.txt") into pr1
    """
    mkdir test
    echo $bin > bin_list.txt
    for b in ${bin}/*.tar.gz
    do
        echo \$b >> bin_list.txt
    done
    bin_ID="${bin.baseName}"
    """
}

process abricate_process{
    publishDir "${params.publish_dir}/test2", mode: 'copy'
    input:
    path bin from bin_dirs_abricate
    output:
    tuple val("bin_ID"), path ("only_bin.txt") into pr2
    """
    echo $bin > only_bin.txt
    bin_ID="${bin.baseName}"
    """
}

workflow {
    

}
