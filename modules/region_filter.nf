#!/bin/bash nextflow 


process REGION_FILTER {
    container = "quay.io/idarahu/finemapping_r:v0.1"

    //cpus 2
    //memory '50 GB'
    //time '20m'
        
    //publishDir "${params.outdir}/${params.prefix2}/standardized", mode: 'copy', pattern: 'standard_*'
    //publishDir "${params.outdir}/${params.prefix2}/$phenotype_id/regions", mode: 'copy', pattern: '*_regions.txt'
    publishDir "${params.outdir}/${params.prefix2}/$phenotype_id/filtered_regions", mode: 'copy', pattern: '*_region.txt'
    publishDir "${params.outdir}/${params.prefix2}/$phenotype_id/filtered_variants", mode: 'copy', pattern: '*_variants.txt'

    input:
    tuple val(chromosome), val(phenotype_id), path(regenie_out), file(bgen_file)
    
    output:
    tuple val(chromosome), val(phenotype_id), path('*_variants.txt'), path('*_region.txt'), file(bgen_file)
    //path 'standard_*'
    //path '*_regions.txt'

    shell:
    '''
    Rscript !{baseDir}/bin/region_filter.R \
     --regenie_out !{regenie_out} \
     --phenotype_id !{phenotype_id} \
     --chr !{chromosome} \
     --p_threshold !{params.p_value} \
     --window_size !{params.window_size}\
     --MHC_start !{params.MHC_start} \
     --MHC_end !{params.MHC_end} \
     --remove_MHC !{params.remove_MHC}
    '''
}
