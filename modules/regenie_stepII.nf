#!/bin/bash nextflow

process REGENIE_STEP_2 {
    container = 'quay.io/eqtlcatalogue/regenie:v3.2.1'
  
    publishDir "${params.outdir}/${params.prefix1}/STEP2/", mode:'copy',\
		pattern: "*.regenie*"
    publishDir "${params.outdir}/${params.prefix1}/logs", mode:'copy',\
		pattern:  "*.log"
     
    
    // Input data    
    input:
    tuple val(chromosome), file(bgen_file), file(bgi_file), file(sample_file)
    path loco_file
    file samples_to_keep
    file pred_list
    file phenotype_file
    val phenotype_list
    file covariate_file   
    
    // Output data
    output:
    tuple val(chromosome), file(bgen_file), file(bgi_file), file(sample_file)
    path "*.regenie*" 
    path "*.log"
    
    shell:         
    '''
    
    nsamples=$(wc -l < !{samples_to_keep}) 
    minMAC=$((nsamples * 2 / 1000))
    
    regenie \
    --step 2 \
    --bgen !{bgen_file} \
    --sample !{sample_file} \
    --ref-first \
    --keep !{samples_to_keep} \
    --phenoFile !{phenotype_file} \
    --phenoColList !{phenotype_list} \
    --covarFile !{covariate_file} \
    --covarColList !{params.covariate_list} \
    --chr !{chromosome} \
    --bsize !{params.bsize} \
    --minMAC $minMAC \
    --minINFO !{params.minINFO} \
    --apply-rint \
    --threads !{task.cpus} \
    --pred !{pred_list} \
    --gz \
    --out chr!{chromosome} \

    '''
}
