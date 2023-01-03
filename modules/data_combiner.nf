#!/bin/bash nextflow

process DATA_COMBINER {
 
    //cpus 4
    //memory '128 GB'
    //time '1h'
        
    publishDir "${params.outdir}/final_data", mode: 'copy', pattern: '*.tsv'
    
    input:
    file coloc5_list
    file coloc3_list
    file clpp_list

    output:
    file 'coloc5_final.tsv'
    file 'coloc3_final.tsv'
    file 'clpp_final.tsv'
    
    shell:
    '''
    head -n 1 !{coloc5_list[0]} > coloc5_final.out
    head -n 1 !{coloc3_list[0]} > coloc3_final.out
    head -n 1 !{clpp_list[0]} > clpp_final.out

    for coloc5_file in !{coloc5_list}
    do
       tail -n+2 -q $coloc5_file >> coloc5_final.out
    done


    for coloc3_file in !{coloc3_list}
    do
       tail -n+2 -q $coloc3_file >> coloc3_final.out
    done


    for clpp_file in !{clpp_list}
    do
       tail -n+2 -q $clpp_file >> clpp_final.out
    done


    cp coloc5_final.out coloc5_final.tsv
    cp coloc3_final.out coloc3_final.tsv
    cp clpp_final.out clpp_final.tsv


    '''
}
