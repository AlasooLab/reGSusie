#!/bin/bash nextflow 


process MANHATTAN_PLOT {
    container = "quay.io/ida_rahu/manhattan_plot_r:latest"
        
    publishDir "${params.outdir}/Manhattan_plots/", mode: 'move', pattern: '*.png'
    publishDir "${params.outdir}/combined_regenie_outputs/", mode: 'copy', pattern: '*.gz'


    input:
    tuple val(phenotype_id), path(regenie_files)
    
    output:
    path('*.gz')
    path('*.png')
    
    shell:
    '''
    
    zcat !{regenie_files[0]} | head -n 1 > !{phenotype_id}_full_regenie.out

    for regenie_file in !{regenie_files} 
    do
      zcat $regenie_file | tail -n +2 >> !{phenotype_id}_full_regenie.out
    done

    Rscript !{baseDir}/bin/Manhattan_plot.R --phenotype_id !{phenotype_id} --regenie_file !{phenotype_id}_full_regenie.out    
    '''
}
