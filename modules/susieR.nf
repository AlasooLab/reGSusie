#!/bin/bash nextflow 


process SUSIER {
    container = "quay.io/idarahu/finemapping_r:v0.1"

    //cpus 2
    //memory '128 GB'
    //time '1h'
        
    publishDir "${params.outdir}/${params.prefix2}/$phenotype_id/susieR", mode: 'copy', pattern: '*.susie.log'
    publishDir "${params.outdir}/${params.prefix2}/coloc5", mode: 'copy', pattern: '*_coloc5.tsv'
    publishDir "${params.outdir}/${params.prefix2}/coloc3", mode: 'copy', pattern: '*_coloc3.tsv'
    publishDir "${params.outdir}/${params.prefix2}/clpp", mode: 'copy', pattern: '*_clpp.tsv'

    input:
    tuple val(phenotype_id), val(region), file(samples), file(filtered_regions), file(ldmatrix)
    
    output:
    file '*_coloc5.tsv'
    file '*_coloc3.tsv'
    file '*_clpp.tsv'
    file '*.susie.log'

    shell:
    '''
    nrows=$(wc -l < !{samples})
    nsamples=$((nrows-2))

    Rscript !{baseDir}/bin/susieR.R \
    --z !{filtered_regions} \
    --ld !{ldmatrix} \
    --log !{phenotype_id}_!{region}.susie.log \
    --pheno !{phenotype_id} \
    --region !{region} \
    --n-samples ${nsamples} \
    --L !{params.max_causal_SNPs} \
    --min-cs-corr 0.5 \
    --n-covariates !{params.n_covariates}
    '''
}
