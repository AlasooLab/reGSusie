#!/bin/bash nextflow

process SEGMENTATION_FAULT_CATCHER {
        
    publishDir "${params.outdir}/segmentation_fault", mode: 'copy', pattern: '*.txt'
    
    input:
    val str

    output:
    path 'regions_list.txt'
    path 'finemapping_again.txt'
    
    shell: 
    '''
    parent_folder=!{baseDir}/!{params.outdir}/!{params.prefix2}
    output_file1="regions_list.txt"
    output_file2="finemapping_again.txt"

    phenotypes=($(find "$parent_folder" -maxdepth 1 -type d ! -name "filtered_regions" -printf "%f\\n" | while read -r phenotype; do
        folder_path="$parent_folder/$phenotype/filtered_regions"
        if [[ -d "$folder_path" ]]; then
            readlink -f "$parent_folder/$phenotype"
        fi
    done))

    for phenotype in "${phenotypes[@]}"; do
        folder_path="$phenotype/filtered_regions"

        for filename in "$folder_path"/"${phenotype##*/}"_*.txt; do
            part=$(echo "$filename" | sed -E 's/.*_([^_]+)_region\\.txt/\\1/')
            z_file="$folder_path/$part.z"

            if [[ "$filename" != *NULL* && ! -f "$z_file" && "$filename" != *finemapping/_*.txt ]]; then
                echo "$filename" >> "$output_file1"

                chr=$(echo "$filename" | sed -E 's/.*_([^_]+)_region\\.txt/\\1/' | cut -d':' -f1)
                regenie_file=!{baseDir}/!{params.outdir}/!{params.prefix1}/"STEP2/chr${chr}_${phenotype##*/}.regenie.gz"
            
                if ! grep -qF "$regenie_file" "$output_file2"; then
                    echo "$regenie_file" >> "$output_file2"
                fi
            fi
        done
    done
    '''
}
