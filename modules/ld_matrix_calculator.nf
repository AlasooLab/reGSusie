#!/bin/bash nextflow

process LD_MATRIX_CALCULATOR {
    container = "quay.io/idarahu/ldmatrix:v0.1"

    //cpus 4
    //memory '128 GB'
    //time '2h'
    //beforeScript 'ulimit -Ss unlimited'
    
    publishDir "${params.outdir}/${params.prefix2}/$phenotype_id/filtered_bgen", mode: 'copy', pattern: '*_filtered*'
    publishDir "${params.outdir}/${params.prefix2}/$phenotype_id/filtered_regions", mode: 'copy', pattern: '*.z'
    publishDir "${params.outdir}/${params.prefix2}/$phenotype_id/ldstore", mode: 'copy', pattern: '*.bcor'
    publishDir "${params.outdir}/${params.prefix2}/$phenotype_id/ldstore", mode: 'copy', pattern: '*.ld*'

    input:
    tuple val(chromosome), val(phenotype_id), val(region), file(variants), file(filtered_regions), file (bgen_file)   
    file sample_file

    output:
    tuple val(phenotype_id), val(region), file(sample_file), file(filtered_regions), file('*.ld*')
    path '*_filtered.bgen'
    path '*_filtered.bgen.bgi'
    path '*.z'
    path '*.bcor'
    
    shell:
    '''
    bgenix -g !{bgen_file} -index
    bgenix -g !{bgen_file} -incl-rsids !{variants} > !{region}_filtered.bgen

    bgenix -g !{region}_filtered.bgen -index
    
    cp !{filtered_regions} !{region}.z

    nrows=$(wc -l < !{sample_file}) 
    nsamples=$((nrows-2))

    echo "z;bgen;bgi;bcor;ld;n_samples" > master.txt
    echo "!{region}.z;!{region}_filtered.bgen;!{region}_filtered.bgen.bgi;!{region}.bcor;!{region}.ld;${nsamples}" >> master.txt   
    
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
