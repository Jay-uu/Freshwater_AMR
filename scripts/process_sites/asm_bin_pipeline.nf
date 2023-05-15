#!/usr/bin/env nextflow
/*
========================================================================================
Multiple co-assemblies of metagenomic soil and water samples
========================================================================================
Trimmed reads have been sorted into different directories for different assembly strategies
----------------------------------------------------------------------------------------
*/

process assemble_water {
    publishDir "${params.publish_dir}/03_megahit", mode: 'copy', pattern: "${sample.baseName}_asm/${sample.baseName}_megahit.fa.gz"
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


process assemble_soil {
    publishDir "${params.publish_dir}/03_megahit", mode: 'copy', pattern: "${sample.baseName}_asm/${sample.baseName}_megahit.fa.gz"
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
    megahit -1 \$read1 -2 \$read2 -o \${sample_ID}_asm --out-prefix \${sample_ID}_megahit --min-contig-len 1000 -t 20 --presets meta-large --min-count 6
    mv \${sample_ID}_asm/\${sample_ID}_megahit.contigs.fa \${sample_ID}_asm/\${sample_ID}_megahit.fa
    gzip \${sample_ID}_asm/\${sample_ID}_megahit.fa
    """
}


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



process binning {
    input:
    tuple val(sample), path(asm), path(mapped_reads)
    output:
    tuple val("$sample"), path("${sample}_bins")
    """
    module load bioinfo-tools MetaBat/2.12.1
    jgi_summarize_bam_contig_depths --outputDepth ${sample}_depth.txt *_sorted.bam
    mkdir ${sample}_bins
    metabat2 -i ${asm} -a ${sample}_depth.txt -o ${sample}_bins/${sample}_bin- --unbinned
    cd ${sample}_bins
    rename bin-. bin- ${sample}_bin-*
    
    #only zero pad if more than 10 bins (if 10 there's an unbinned too)
    bin_count=`ls *.fa | wc -l`
    if ((\$bin_count > 10));
    then
        char_count=`echo -n \$bin_count | wc -c`
        v=9
        for ((i=1; i < \$char_count-1; i++)); do v+=9; done
        for (( i=1; i <= \$v; i++ )); do mv -- ${sample}_bin-\${i}.fa ${sample}_bin-\$(printf '%0*d\n' \$char_count \$i).fa; done
    fi
    cd ..
    """
}

process rename_contigs {
    publishDir "${params.publish_dir}/05_metabat", mode: 'copy'
    conda "${params.envs_dir}/biopython-1.78"
    input:
    tuple val(sample), path(bins_dir)
    output:
    tuple val("$sample"), path("$bins_dir")
    """
    #!/usr/bin/env python
    from Bio import SeqIO
    import math
    import os
    for bin in os.listdir("$bins_dir"):
        all_contigs = list(SeqIO.parse("${bins_dir}"+'/'+bin, "fasta"))
        num_len = int(math.log10(len(all_contigs))+1)
        bin_name = bin.split('.',1)[0]
        for i,s in enumerate(all_contigs):
            s.id =f"{bin_name}_{i+1:0{num_len}}" #renames contigs
        with open("${bins_dir}"+'/'+bin, "w") as output_handle:
            SeqIO.write(all_contigs, output_handle, "fasta")
    """
}



process prokka_annot {
    publishDir "${params.publish_dir}/06_prokka/", mode: 'copy'
    input:
    tuple val(sample), path(bins_dir)
    output:
    tuple val("$sample"), path("${sample}_prokka")
    """
    module load bioinfo-tools prokka/1.45-5b58020
    mkdir ${sample}_prokka
    for bin in $bins_dir/*
    do
    bin_name=`basename \$bin .fa`
    prokka --outdir \${bin_name} --prefix \$bin_name \$bin --cpus 0
    gzip \${bin_name}/\$bin_name.*
    tar -czf \${bin_name}.tar.gz \$bin_name
    mv \${bin_name}.tar.gz ${sample}_prokka/
    done
    """
}

process rename_prokka_out {
    publishDir "${params.publish_dir}/07_bins/", mode: 'copy'
    conda "${params.envs_dir}/biopython-1.78"
    input:
    tuple val(sample), path(prokka_dir)
    output:
    tuple val("$sample"), path("$sample")
    """
    #want fna.gz, ffn.gz, faa.gz, gff.gz
    mkdir $sample
    for bin in $prokka_dir/*.tar.gz
    do
    bin_name=`basename \$bin .tar.gz`
    mkdir ${sample}/\$bin_name
    tar -xf \$bin -C ${sample}/\$bin_name --wildcards '*.fna.gz' --strip-components=1
    tar -xf \$bin -C ${sample}/\$bin_name --wildcards '*.ffn.gz' --strip-components=1
    tar -xf \$bin -C ${sample}/\$bin_name --wildcards '*.faa.gz' --strip-components=1
    tar -xf \$bin -C ${sample}/\$bin_name --wildcards '*.gff.gz' --strip-components=1
    done
    #only ffn and faa are going to need to be renamed
    for bin in ${sample}/*
    do
    cd \$bin
    gunzip *.f*.gz
    python /crex/proj/fume/nobackup/private/jay/Freshwater_AMR/scripts/process_sites/rename_contigs.py -i *.ffn -o *.ffn
    python /crex/proj/fume/nobackup/private/jay/Freshwater_AMR/scripts/process_sites/rename_contigs.py -i *.faa -o *.faa
    gzip *.f*
    cd ../..
    tar -czf \${bin}.tar.gz \$bin
    rm -r \$bin
    done
    """
}




workflow {
    params.publish_dir="/proj/fume/private/jay/process_sites"
    params.envs_dir="/proj/fume/nobackup/private/jay/Freshwater_AMR/conda_envs"
    
    //groups are directories with multiple trimmed PE reads, while singles are trimmed PE reads for individual samples
    params.water_groups="/crex/proj/fume/private/jay/process_sites/02_fastp/group_sets/water"
    params.soil_groups="/crex/proj/fume/private/jay/process_sites/02_fastp/group_sets/soil"
    params.water_singles="/crex/proj/fume/private/jay/process_sites/02_fastp/water"
    params.soil_singles="/crex/proj/fume/private/jay/process_sites/02_fastp/soil"
    
    //for testing
    //water_group = Channel.fromPath( "${params.water_groups}/Abisko", type:'dir' )
    //water = Channel.fromPath( "${params.water_singles}/Sample-A1/", type:'dir'  )
    //assemble_water(water) | map_reads_asm | binning | rename_contigs | prokka_annot

    water_coasm = Channel.fromPath( "${params.water_groups}/*", type:'dir' )
    soil_coasm = Channel.fromPath( "${params.soil_groups}/*", type:'dir' )
    water_single = Channel.fromPath( "${params.water_singles}/*", type:'dir'  )
    soil_single = Channel.fromPath( "${params.soil_singles}/*", type:'dir'  )
    
    //concatenate channels
    water = water_coasm.concat(water_single)
    soil = soil_coasm.concat(soil_single)
    assemble_water(water)
    assemble_soil(soil)
    all_samps = assemble_water.out.concat(assemble_soil.out)
    map_reads_asm(all_samps) | binning | rename_contigs | prokka_annot | rename_prokka_out

//assemble_water(water) | map_reads_asm | binning //| rename_contigs | prokka_annot
  //  assemble_soil(soil_coasm) | map_reads_asm | binning //| rename_contigs | prokka_annot
    

}

