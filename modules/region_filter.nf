#!/bin/bash nextflow 


process REGION_FILTER {
    container = "quay.io/idarahu/finemapping_r:v0.1"
        
    publishDir "${params.outdir}/${params.prefix2}/$phenotype_id/filtered_regions", mode: 'copy', pattern: '*_region.txt'


    input:
    tuple val(chromosome), file(bgen_file), file(bgi_file), file(sample_file), val(phenotype_id), path(regenie_out)
    
    output:
    tuple val(chromosome), val(phenotype_id), file(bgen_file), file(bgi_file), file(sample_file), path('*_region.txt')
    
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
