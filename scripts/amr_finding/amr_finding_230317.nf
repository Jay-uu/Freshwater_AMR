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

process checkm_process{
    publishDir "${params.publish_dir}/01_checkm", mode: 'copy'
    input:
    path sample
    output:
    tuple val("${sample.baseName}"), path("${sample.baseName}_summary.txt")
    """
    echo '===== Start ${sample.baseName} ====='
    date
    mkdir tmp_faa
    for arch in $sample/*.tar.gz
    do
    tar -xf \$arch -C tmp_faa --wildcards '*.faa.gz' --strip-components=3
    done
    gunzip tmp_faa/*.faa.gz
    checkm lineage_wf --genes -t 20 -x faa tmp_faa ${sample.baseName}_res
    checkm qa ${sample.baseName}_res/lineage.ms ${sample.baseName}_res > ${sample.baseName}_summary.txt
    rm -r tmp_faa
    date
    echo '===== Done ${sample.baseName} ====='
    """

}

process parse_checkm_summary{
    publishDir  "${params.publish_dir}/01_checkm", mode: 'copy'
    input:
    tuple val(sample_basename), path(summary)
    output:
    tuple val("${sample_basename}"), path("${sample_basename}_summary.json")
    """
    python3 /home/jay/master_thesis/scripts/parse_checkm_output_edit.py $summary > ${sample_basename}_summary.json
    """

}

process run_abricate{
    publishDir "${params.publish_dir}/02_abricate", mode: 'copy'
    input:
    path sample
    output:
    tuple val("${sample.baseName}"), path("${sample.baseName}_abricate.tab")
    """
    echo '===== Start ${sample.baseName} ====='
    date
    mkdir ${sample.baseName}_bins
    for bin in $sample/*.tar.gz
    do
    tar -xf \$bin -C ${sample.baseName}_bins --wildcards '*.ffn.gz' --strip-components=3
    done
    abricate --db card --minid 50 ${sample.baseName}_bins/*.ffn.gz > ${sample.baseName}_abricate.tab
    date
    echo '===== Done ${sample.baseName} ====='
    """
}

process run_rgi{
    publishDir "${params.publish_dir}/03_rgi", mode: 'copy'
    input:
    path sample
    output:
    tuple val("${sample.baseName}"), path("${sample.baseName}_res")
    """
    echo '===== Start ${sample.baseName} ====='
    date
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
    echo "running bin: \$file"
    BIN_ID=`basename \$file .faa`
    rgi main --input_sequence \$file --output_file \${sample_ID}_res/\${BIN_ID} -t protein -n 23 -a DIAMOND --low_quality --clean
    done
    date
    echo '===== Done ${sample.baseName} ====='
    """
}

process run_gtdbtk{
    publishDir "${params.publish_dir}/04_gtdbtk", mode: 'copy'
    input:
    path sample
    output:
    tuple val("${sample.baseName}"), path("${sample.baseName}_result")
    """
    echo '===== Running ${sample.baseName} ====='
    date
    mkdir ${sample.baseName}_bins
    mkdir ${sample.baseName}_res
    for bin in $sample/*.tar.gz
    do
    tar -xf \$bin -C ${sample.baseName}_bins --wildcards '*.fna.gz' --strip-components=3
    done
    gtdbtk classify_wf --genome_dir ${sample.baseName}_bins --out_dir ${sample.baseName}_result --cpus 20 -x gz
    date
    echo '===== Done ${sample.baseName} ====='
    """

}

process combine_all{
    publishDir "${params.publish_dir}/05_combined_amr", mode: 'copy'
    input:
    tuple val(sample), path(checkm_res), path(abr_res), path(rgi_dir), path(gtdbtk_dir)
    output:
    path "${sample}.json"
    """
    python3 /home/jay/master_thesis/scripts/parse_amr_res_nf.py -rgi $rgi_dir -abr $abr_res -chm $checkm_res -tax ${gtdbtk_dir}/gtdbtk.bac120.summary.tsv -o ${sample}.json       
    """

}
workflow {

    Channel
        //.fromPath("/home/moritz/data/data_submit/bins/ErkenS*", type:'dir')
        .fromPath("/home/jay/data/test_bins/*", type: 'dir')
        .multiMap { bin_dirs -> bin_dirs_checkm: bin_dirs_abricate: bin_dirs_rgi: bin_dirs_gtdbtk: bin_dirs }
        .set { dirs }
    //results directory
    params.publish_dir='/home/jay/data/test_amr_finding/'

    checkm_process(dirs.bin_dirs_checkm) | parse_checkm_summary
    run_abricate(dirs.bin_dirs_abricate)
    run_rgi(dirs.bin_dirs_rgi)
    run_gtdbtk(dirs.bin_dirs_gtdbtk)

    parse_checkm_summary.out \
        | combine( run_abricate.out, by: 0) \
        | combine( run_rgi.out, by: 0) \
        | combine( run_gtdbtk.out, by: 0) \
        | combine_all \
}
