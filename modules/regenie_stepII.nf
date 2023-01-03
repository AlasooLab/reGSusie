#!/bin/bash nextflow
//nextflow.enable.dsl=2

process REGENIE_STEP_2 {
    container = 'quay.io/eqtlcatalogue/regenie:v3.2.1'

    // Process directives
    //cpus 16
    //memory '4 GB'
    //time '1h'
    

    publishDir "${params.outdir}/${params.prefix1}/STEP2/${phenotype_id}/", mode:'copy',\
		pattern: "*.regenie*"
    publishDir "${params.outdir}/${params.prefix1}/logs", mode:'copy',\
		pattern:  "*.log"
     
    
    // Input data    
    input:
    tuple val(chromosome), file(bgen_file), val(phenotype_id), path(loco_file)
    file sample_file
    file pred_list
    file phenotype_file
    file covariate_file   
    
    // Output data
    output:
    tuple val(chromosome), val(phenotype_id), path("*.regenie*"), file(bgen_file)
    path "*.log"
    
    shell:         
    '''
    regenie \
    --step 2 \
    --bgen !{bgen_file} \
    --sample !{sample_file} \
    --ref-first \
    --phenoFile !{phenotype_file} \
    --phenoColList !{phenotype_id} \
    --covarFile !{covariate_file} \
    --covarColList !{params.covariate_list} \
    --chr !{chromosome} \
    --bsize 1000 \
    --apply-rint \
    --threads !{task.cpus} \
    --pred !{pred_list} \
    --gz \
    --out chr!{chromosome} \

    '''
}
