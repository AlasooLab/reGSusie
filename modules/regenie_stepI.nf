#!/bin/bash nextflow
//nextflow.enable.dsl=2

process REGENIE_STEP_1 {
    container = 'quay.io/eqtlcatalogue/regenie:v3.2.1'
       
    publishDir "${params.outdir}/${params.prefix1}/STEP1", mode:'copy',\
		pattern: "${params.prefix1}*"
    publishDir "${params.outdir}/${params.prefix1}/logs", mode:'copy',\
		pattern:  ".command.out", saveAs: { filename -> "$params.prefix1-step1.out" }
    publishDir "${params.outdir}/${params.prefix1}/STEP1", mode:'copy',\
		pattern:  ".command.sh", saveAs: { filename -> "$params.prefix1-step1.sh" }

    // Input data
    input:
    file bgen_file
    file sample_file
    file phenotype_file
    val phenotype_list
    file covariate_file

    // Output data
    output:
    path "${params.prefix1}_pred.list"
    file "${params.prefix1}*.loco.gz"
    path ".command.out"
    path ".command.sh"

    shell:
    '''
    
    regenie \
    --step 1 \
    --bgen !{bgen_file} \
    --sample !{sample_file} \
    --ref-first \
    --phenoFile !{phenotype_file} \
    --phenoColList !{phenotype_list} \
    --covarFile !{covariate_file} \
    --covarColList !{params.covariate_list} \
    --use-relative-path \
    --bsize !{params.bsize} \
    --lowmem \
    --apply-rint \
    --lowmem-prefix !{params.prefix1}_tmp_rg \
    --gz \
    --threads !{task.cpus} \
    --out !{params.prefix1} \
    
    '''
}
