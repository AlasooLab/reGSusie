#!/bin/bash nextflow
nextflow.enable.dsl=2

// Importing modules
include { REGENIE_STEP_1 } from './modules/regenie_stepI'
include { REGENIE_STEP_2 } from './modules/regenie_stepII'
include {REGION_FILTER} from './modules/region_filter'
include {LD_MATRIX_CALCULATOR} from './modules/ld_matrix_calculator'
include {SUSIER} from './modules/susieR'
include {DATA_COMBINER} from './modules/data_combiner'

// Workflow
workflow {
    bgen_step1_ch = Channel.fromPath(params.step1_bgen).collect()
    sample_step1_ch = Channel.fromPath(params.step1_sample).collect()

    phenotype_ch = Channel.fromPath(params.phenotype_file).collect()
    phenotype_list_ch = Channel.value(file(params.phenotype_list).text)
    covariate_ch = Channel.fromPath(params.covariate_file).collect()
    
    
    REGENIE_STEP_1(bgen_step1_ch, sample_step1_ch,  phenotype_ch, phenotype_list_ch, covariate_ch)
    
    
    bgen_step2_ch = Channel.fromPath(params.step2_bgen).flatten()
						       .map { file ->
                                                              tuple(file.simpleName.split('chr')[1].split('_')[0], file)
                                                            }
    bgi_step2_ch = Channel.fromPath(params.step2_bgi).flatten()
						       .map { file ->
                                                              tuple(file.simpleName.split('chr')[1].split('_')[0], file)
                                                            }
    sample_step2_ch = Channel.fromPath(params.step2_sample).flatten()
                      .map { file ->
                                                              tuple(file.simpleName.split('chr')[1].split('_')[0], file)
                                                            }
    
    bgen_sample_ch = bgen_step2_ch.join(bgi_step2_ch, by: 0).join(sample_step2_ch, by: 0)
    samples_to_keep_ch = Channel.fromPath(params.samples_to_keep).collect()
    
    pred_list_ch = REGENIE_STEP_1.out[0]
    loco_file_ch = REGENIE_STEP_1.out[1]
    
    
    REGENIE_STEP_2(bgen_sample_ch, loco_file_ch, samples_to_keep_ch, pred_list_ch, phenotype_ch, phenotype_list_ch, covariate_ch)
    
    
    chr_bgen_sample_ch = REGENIE_STEP_2.out[0]
    regenie_file_ch = REGENIE_STEP_2.out[1].flatten()
                      .map { file -> tuple(file.simpleName.split('_')[1..-1].join('_'), file)
                      }

    region_ch = chr_bgen_sample_ch.combine(regenie_file_ch)
    
    
    REGION_FILTER(region_ch)
    

    ld_input_ch = REGION_FILTER.out[0].transpose()
                  .filter { it[5].name.split('_')[-2] != 'NULL' }
    				      .map { chr, pheno_id, bgen_file, bgi_file, sample_file, region_file ->
                         tuple(chr, pheno_id, region_file.name.split('_')[-2], region_file, bgen_file, bgi_file, sample_file) 
    				           }

    samples_ld_incl_ch = Channel.fromPath(params.samples_ld_incl).collect()                                  
    
    
    LD_MATRIX_CALCULATOR(ld_input_ch, samples_ld_incl_ch)
    
    
    SUSIER(LD_MATRIX_CALCULATOR.out[0])
    
    
    clpp_ch = SUSIER.out[2].filter {it.name.split('_')[-2] != 'NULL'}.collect()
    
    
    DATA_COMBINER(SUSIER.out[0].collect(), SUSIER.out[1].collect(), clpp_ch) 
}
