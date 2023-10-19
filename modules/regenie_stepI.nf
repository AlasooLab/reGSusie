#!/bin/bash nextflow

mode = params.step1_input_format
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
  file phenotype_file
  val phenotype_list
  file covariate_file

  // Output data
  output:
  path "${params.prefix1}_pred.list"
  file "${params.prefix1}*.loco.gz"
  path ".command.out"
  path ".command.sh"

  script:
  
  if ( mode == 'bed' ) 
    """
    regenie \
    --step 1 \
    --bed ${params.step1_bed} \
    --ref-first \
    --phenoFile ${phenotype_file} \
    --phenoColList ${phenotype_list} \
    --covarFile ${covariate_file} \
    --covarColList ${params.covariate_list} \
    --use-relative-path \
    --bsize ${params.bsize} \
    --lowmem \
    --apply-rint \
    --lowmem-prefix ${params.prefix1}_tmp_rg \
    --gz \
    --threads ${task.cpus} \
    --out ${params.prefix1} \
    
    """
    
  else if ( mode == 'pgen' ) 
    """
    regenie \
    --step 1 \
    --pgen ${params.step1_pgen} \
    --ref-first \
    --phenoFile ${phenotype_file} \
    --phenoColList ${phenotype_list} \
    --covarFile ${covariate_file} \
    --covarColList ${params.covariate_list} \
    --use-relative-path \
    --bsize ${params.bsize} \
    --lowmem \
    --apply-rint \
    --lowmem-prefix ${params.prefix1}_tmp_rg \
    --gz \
    --threads ${task.cpus} \
    --out ${params.prefix1} \
    
    """
    
  else if ( mode == 'bgen' ) 
    """
    regenie \
    --step 1 \
    --bgen ${params.step1_bgen} \
    --sample ${params.step1_sample} \
    --ref-first \
    --phenoFile ${phenotype_file} \
    --phenoColList ${phenotype_list} \
    --covarFile ${covariate_file} \
    --covarColList ${params.covariate_list} \
    --use-relative-path \
    --bsize ${params.bsize} \
    --lowmem \
    --apply-rint \
    --lowmem-prefix ${params.prefix1}_tmp_rg \
    --gz \
    --threads ${task.cpus} \
    --out ${params.prefix1} \
    
    """
  
  else
    error "Invalid input format: ${params.step1_input_format}"
    
}