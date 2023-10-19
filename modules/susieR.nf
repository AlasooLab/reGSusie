#!/bin/bash nextflow 


process SUSIER {
    container = "quay.io/idarahu/finemapping_r:v0.1"
        
    publishDir "${params.outdir}/${params.prefix2}/$phenotype_id/susieR", mode: 'copy', pattern: '*.susie.log'
    publishDir "${params.outdir}/susieR/coloc5", mode: 'copy', pattern: '*_coloc5.tsv'
    publishDir "${params.outdir}/susieR/coloc3", mode: 'copy', pattern: '*_coloc3.tsv'
    publishDir "${params.outdir}/susieR/clpp", mode: 'copy', pattern: '*_clpp.tsv'

    input:
    tuple val(phenotype_id), val(region), file(samples), file(filtered_regions), file(ldmatrix)
    
    output:
    tuple val(phenotype_id), path('*_coloc5.tsv')
    tuple val(phenotype_id), path('*_coloc3.tsv')
    tuple val(phenotype_id), path('*_clpp.tsv')
    path '*.susie.log'

    shell:
    '''
    nsamples=$(wc -l < !{samples})

    Rscript !{baseDir}/bin/susieR.R \
    --z !{filtered_regions} \
    --ld !{ldmatrix} \
    --log !{phenotype_id}_!{region}.susie.log \
    --pheno !{phenotype_id} \
    --region !{region} \
    --n-samples ${nsamples} \
    --L !{params.max_causal_SNPs} \
    --min-cs-corr 0.5 \
    --n-covariates !{params.n_covariates} \
    --GRCh !{params.GRCh}
    '''
}
