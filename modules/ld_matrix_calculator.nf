#!/bin/bash nextflow

process LD_MATRIX_CALCULATOR {
    container = "quay.io/idarahu/ldmatrix:v0.1"
    
    publishDir "${params.outdir}/${params.prefix2}/$phenotype_id/filtered_regions", mode: 'copy', pattern: '*.z'
    publishDir "${params.outdir}/${params.prefix2}/$phenotype_id/ldstore", mode: 'copy', pattern: '*.bcor'
    publishDir "${params.outdir}/${params.prefix2}/$phenotype_id/ldstore", mode: 'copy', pattern: '*.ld*'

    input:
    tuple val(chromosome), val(phenotype_id), val(region), file(filtered_regions), file(bgen_file), file(bgi_file), file(sample_file)     
    file sample_file_incl

    output:
    tuple val(phenotype_id), val(region), file(sample_file_incl), file(filtered_regions), file('*.ld*')
    path '*.z'
    path '*.bcor'
    
    shell:
    '''  
    cp !{filtered_regions} !{region}.z

    nsamples=$(wc -l < !{sample_file_incl}) 

    echo "z;bgen;bgi;bcor;ld;n_samples;sample;incl" > master.txt
    echo "!{region}.z;!{bgen_file};!{bgi_file};!{region}.bcor;!{region}.ld;${nsamples};!{sample_file};!{sample_file_incl}" >> master.txt   
    
    ldstore \
    --in-files master.txt \
    --write-bcor \
    --read-only-bgen \
    --n-threads !{task.cpus}

    ldstore \
    --in-files master.txt \
    --bcor-to-text

    gzip !{region}.ld
    '''
}
