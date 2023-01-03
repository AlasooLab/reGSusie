#!/bin/bash nextflow 


process MANHATTAN_PLOT {
    container = "quay.io/eqtlcatalogue/regenie_manhattan:v22.09.1"

    cpus 2
    memory '12 GB'
    time '10m'
        
    publishDir "${params.outdir}/${params.prefix}/Manhattan_plots", mode: 'move', pattern: "*.png"
    publishDir "${params.outdir}/${params.prefix}/logs", mode:'copy',\
                pattern: ".command.out", saveAs: { filename -> "$phenotype_id-Manhattan_plots.out" }
    
    input:
    tuple val(phenotype_id), path(regenie_out)

    output:
    path "*.png"
    path ".command.out"

    shell:
    '''
    Rscript !{baseDir}/bin/Manhattan_plot.R \
     --file !{regenie_out} \
     --phenotype_id !{phenotype_id} \
     --out !{phenotype_id} 
    '''
}

