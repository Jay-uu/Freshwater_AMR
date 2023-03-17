#!/usr/bin/env nextflow
/*
========================================================================================
AMR-finding in multiple metagenomic assembled lake samples
========================================================================================
With quality eval and taxonomic identification
----------------------------------------------------------------------------------------
*/


// Get archives with bins
//bin_dirs=Channel.fromPath("/home/moritz/data/data_submit/bins/*", type:'dir')
//start with one bin. I'm thinking ErkenSummer because then I can compare and see if the results are the same!
bin_dirs=Channel.fromPath("/home/moritz/data/data_submit/bins/ErkenSum*", type:'dir')
bin_dirs.into { bin_dirs_checkm; bin_dirs_abricate; bin_dirs_rgi; bin_dirs_gtdbtk }

params.publish_dir='/home/jay/data/amr_finding'

process checkm_process{
    publishDir "${params.publish_dir}/01_checkm", mode: 'copy'
    input:
    path bin from bin_dirs_checkm
    output:
    tuple val("${bin.baseName}"), path("${bin.baseName}_summary.txt") into checkm_summary
    """
    mkdir tmp_faa
    for arch in $bin/*.tar.gz
    do
    tar -xf \$arch -C tmp_faa --wildcards '*.faa.gz' --strip-components=3
    done
    gunzip tmp_faa/*.faa.gz
    checkm lineage_wf --genes -t 20 -x faa tmp_faa ${bin.baseName}_res
    checkm qa ${bin.baseName}_res/lineage.ms ${bin.baseName}_res > ${bin.baseName}_summary.txt
    rm -r tmp_faa
    """

}

process parse_checkm_summary{
    publishDir  "${params.publish_dir}/01_checkm", mode: 'copy'
    input:
    tuple val(bin_basename), path(summary) from checkm_summary
    output:
    tuple val("${bin_basename}"), path("${bin_basename}_summary.json")
    """
    python3 /home/jay/master_thesis/scripts/parse_checkm_output_edit.py $summary > ${bin_basename}_summary.json
    """

}

process run_abricate{
    publishDir "${params.publish_dir}/02_abricate", mode: 'copy'
    input:
    path sample from bin_dirs_abricate
    output:
    tuple val("${sample.baseName}"), path("${sample.baseName}_abricate.tab") into abricate_res
    """
    mkdir ${sample.baseName}_bins
    for bin in $sample/*.tar.gz
    do
    tar -xf \$bin -C ${sample.baseName}_bins --wildcards '*.ffn.gz' --strip-components=3
    done
    abricate --db card --minid 50 ${sample.baseName}_bins/*.ffn.gz > ${sample.baseName}_abricate.tab
    """
}

process run_rgi{
    publishDir "${params.publish_dir}/03_rgi", mode: 'copy'
    input:
    path sample from bin_dirs_rgi
    output:
    tuple val("${sample.baseName}"), path("${sample_ID}_res") into rgi_res
    """
    sample_ID="${sample.baseName}"
    mkdir \${sample_ID}_bins
    for bin in $sample/*.tar.gz
    do
    tar -xf \$bin -C \${sample_ID}_bins --wildcards '*.faa.gz' --strip-components=3
    done
    
    gunzip \${sample_ID}_bins/*.faa.gz
    sed -i 's/[*]//g' \${sample_ID}_bins/*.faa
    mkdir \${sample_ID}_res
    for file in \${sample_ID}_bins/*
    do
    BIN_ID=`basename \$file .faa`
    rgi main --input_sequence \$file --output_file \${sample_ID}_res/\${BIN_ID} -t protein -n 23 -a DIAMOND --low_quality --clean
    done
    """
}

process run_gtdbtk{
    publishDir "${params.publish_dir}/04_gtdbtk", mode: 'copy'
    input:
    path sample from bin_dirs_gtdbtk
    output:
    tuple val("${sample.baseName}"), path("${sample.baseName}_res") into gtdbtk_res
    """
    mkdir ${sample.baseName}_bins
    mkdir ${sample.baseName}_res
    for bin in $sample/*.tar.gz
    do
    tar -xf \$bin -C ${sample.baseName}_bins --wildcards '*.fna.gz' --strip-components=3
    done
    gtdbtk classify_wf --genome_dir ${sample.baseName}_bins --out_dir ${sample.baseName}_res --cpus 20 -x gz
    """

}
