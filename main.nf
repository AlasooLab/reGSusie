#!/bin/bash nextflow
nextflow.enable.dsl=2

// Importing modules
include { REGENIE_STEP_1 } from './modules/regenie_stepI'
include { REGENIE_STEP_2 } from './modules/regenie_stepII'
//include { MANHATTAN_PLOT } from './modules/manhattan_plot'
include {REGION_FILTER} from './modules/region_filter'
include {LD_MATRIX_CALCULATOR} from './modules/ld_matrix_calculator'
include {SUSIER} from './modules/susieR'
include {DATA_COMBINER} from './modules/data_combiner'

// Workflow
workflow {
    bgen_step1_ch = Channel.fromPath(params.step1_bgen).collect()
    sample_ch = Channel.fromPath(params.sample).collect()
    phenotype_ch = Channel.fromPath(params.phenotype_file).collect()
    phenotype_list_ch = Channel.value(file(params.phenotype_list).text)
    covariate_ch = Channel.fromPath(params.covariate_file).collect()
    
    REGENIE_STEP_1(bgen_step1_ch, sample_ch,  phenotype_ch, phenotype_list_ch, covariate_ch)
    
    bgen_step2_ch = Channel.fromPath(params.step2_bgen).flatten()
						       .map { file ->
                                                              tuple(file.simpleName.split('chr')[1], file)
                                                            }
    
    pred_list_ch = REGENIE_STEP_1.out[0]
    pheno_id_ch = pred_list_ch.splitCsv(sep: ' ')
    loco_file_ch = REGENIE_STEP_1.out[1].flatten().map { loco_file -> 
                                                         tuple(loco_file, loco_file.name).collect()
                                                       }
    pheno_loco_ch = pheno_id_ch.join(loco_file_ch, by: 1).map { name, pheno_id, loco_file ->
                                                                tuple(pheno_id, loco_file)
                                                              }
    step2_input_ch = bgen_step2_ch.combine(pheno_loco_ch)
    REGENIE_STEP_2(step2_input_ch, sample_ch,  pred_list_ch,  phenotype_ch, covariate_ch)
   
    //MANHATTAN_PLOT(REGENIE_STEP_2.out[0])
    
    REGION_FILTER(REGENIE_STEP_2.out[0])
    
    //control_ch = REGION_FILTER.out[0].transpose() //map {chr, pheno_id, variants_file, region_file, bgen_file -> variants_file.name}
    

    //control_ch.view()


    ld_input_ch = REGION_FILTER.out[0].transpose()
                                      .filter { it[2].name.split('_')[-2] != 'NULL' }
    				      .map { chr, pheno_id, variants_file, region_file, bgen_file ->
                                             tuple(chr, pheno_id, variants_file.name.split('_')[-2], variants_file, region_file, bgen_file) 
    				           }
                                      //.join(bgen_step2_ch, by: 0)

    //ld_input_ch.view() 
    LD_MATRIX_CALCULATOR(ld_input_ch, sample_ch)

    SUSIER(LD_MATRIX_CALCULATOR.out[0])
    
    clpp_ch = SUSIER.out[2].filter {it.name.split('_')[-2] != 'NULL'}.collect()
    //clpp_ch.view()
    DATA_COMBINER(SUSIER.out[0].collect(), SUSIER.out[1].collect(), clpp_ch)   
}
