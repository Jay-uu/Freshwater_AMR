#!/usr/bin/env nextflow
/*
========================================================================================
Combining bin quality, AMR gene hits, taxonomic identity etc. for contigs with hits 
========================================================================================

----------------------------------------------------------------------------------------
*/


process combine_all{
    publishDir "${params.input_dir}/06_combined_contig_amr", mode: 'copy'
    conda "${params.envs_dir}/pandas-1.5.3"
    input:
    tuple val(sample), path(checkm_res), path(abr_res), path(rgi_dir), path(gtdbtk_dir)
    output:
    path "${sample}.json"
    path "${sample}_aro_info.csv"
    """
    python3 /home/jay/master_thesis/scripts/parse_amr_contigs.py -rgi $rgi_dir -abr $abr_res -chm $checkm_res -tax ${gtdbtk_dir}/gtdbtk.bac120.summary.tsv -o ${sample}.json -m ${sample}_aro_info.csv       
    """

}

workflow {

    params.input_dir = "/proj/fume/private/jay/amr_finding"
    params.envs_dir = "/proj/fume/nobackup/private/jay/Freshwater_AMR/conda_envs"
    rem = ~/_.*/
    checkm = Channel
                .fromPath("${params.input_dir}/01_checkm/*.json", type:'file')
                .map { file -> tuple((file.baseName - rem), file)  }
    abricate = Channel
                .fromPath("${params.input_dir}/02_abricate/*", type:'file')
                .map { file -> tuple((file.baseName - rem), file)  }
    rgi = Channel
                .fromPath("${params.input_dir}/03_rgi/*", type:'dir')
                .map { dir -> tuple((dir.baseName - rem), dir)  }
    gtdbtk = Channel
                .fromPath("${params.input_dir}/04_gtdbtk/*", type:'dir')
                .map { dir -> tuple((dir.baseName - rem), dir)  }
    checkm \
        | combine( abricate, by: 0) \
        | combine( rgi, by: 0) \
        | combine( gtdbtk, by: 0) \
        | combine_all \
}
