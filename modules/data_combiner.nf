#!/bin/bash nextflow

process DATA_COMBINER {
        
    publishDir "${params.outdir}/combined_SuSiE_outputs/${file_type}_final/", mode: 'copy', pattern: '*.tsv.gz'
    
    input:
    tuple val(phenotype_id), path(susie_files), val(file_type)

    output:
    file '*.tsv.gz'
    
    shell:
    '''
    
    cat !{susie_files[0]} | head -n 1 > header.txt

    chromosome_column=$(awk -v FS="\t" 'NR==1{for (i=1; i<=NF; i++) if ($i == "chromosome") print i}' header.txt)
    position_column=$(awk -v FS="\t" 'NR==1{for (i=1; i<=NF; i++) if ($i == "position") print i}' header.txt)

    cat !{susie_files[0]} | head -n 1 > !{phenotype_id}_!{file_type}_final.out

    for susie_file in !{susie_files} 
    do
      cat $susie_file | tail -n +2 >> !{phenotype_id}_!{file_type}_final.out
    done
    sort -t$'\t' -k$chromosome_column,$chromosome_column -k$position_column,$position_column -o !{phenotype_id}_!{file_type}_final.out !{phenotype_id}_!{file_type}_final.out

    gzip -c !{phenotype_id}_!{file_type}_final.out > !{phenotype_id}_!{file_type}_final.tsv.gz

    '''
}
