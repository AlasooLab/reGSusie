#!/bin/bash nextflow
nextflow.enable.dsl=2

// Importing modules
include { REGENIE_STEP_1 } from './modules/regenie_stepI'
include { REGENIE_STEP_2 } from './modules/regenie_stepII'
include { MANHATTAN_PLOT } from './modules/Manhattan_plot'
include { REGION_FILTER } from './modules/region_filter'
include { LD_MATRIX_CALCULATOR } from './modules/ld_matrix_calculator'
include { SEGMENTATION_FAULT_CATCHER } from './modules/segmentation_fault_catcher'
include { SUSIER } from './modules/susieR'
include { DATA_COMBINER as DATA_COMBINER1 } from './modules/data_combiner'
include { DATA_COMBINER as DATA_COMBINER2 } from './modules/data_combiner'
include { DATA_COMBINER as DATA_COMBINER3 } from './modules/data_combiner'

// Workflow
workflow {
    
  bgen_step2_ch = Channel.fromPath(params.step2_bgen).flatten()
                         .map { file -> tuple(file.simpleName.split('chr')[1].split('_')[0], file) }
                         
  bgi_step2_ch = Channel.fromPath(params.step2_bgi).flatten()
                        .map { file -> tuple(file.simpleName.split('chr')[1].split('_')[0], file) }
                        
  sample_step2_ch = Channel.fromPath(params.step2_sample).flatten()
                           .map { file -> tuple(file.simpleName.split('chr')[1].split('_')[0], file) }
  
  // This channel is also necessary when running solely the LDstore and SuSiE components of the pipeline.                           
  bgen_sample_ch = bgen_step2_ch.join(bgi_step2_ch, by: 0).join(sample_step2_ch, by: 0)
  
  if ( params.full_pipeline ) {
  
    phenotype_ch = Channel.fromPath(params.phenotype_file).collect()
    phenotype_list_ch = Channel.value(file(params.phenotype_list).text)
    covariate_ch = Channel.fromPath(params.covariate_file).collect()
    
      
    REGENIE_STEP_1(phenotype_ch, phenotype_list_ch, covariate_ch)
    
    
    samples_to_keep_ch = Channel.fromPath(params.samples_to_keep).collect()
    
    pred_list_ch = REGENIE_STEP_1.out[0]
    loco_file_ch = REGENIE_STEP_1.out[1]
    
    
    REGENIE_STEP_2(bgen_sample_ch, loco_file_ch, samples_to_keep_ch, pred_list_ch, phenotype_ch, phenotype_list_ch, covariate_ch)
    
      
    chr_bgen_sample_ch = REGENIE_STEP_2.out[0]
    regenie_file_ch = REGENIE_STEP_2.out[1].flatten()
                      .map { file -> tuple(file.simpleName.split('_')[1..-1].join('_'), file)
                           }
    
    region_ch = chr_bgen_sample_ch.combine(regenie_file_ch)
                                  .filter { it[0] == it[-1].simpleName.split('chr')[1].split('_')[0] }
    
    grouped_regenie_ch = regenie_file_ch.groupTuple()
    
    
    MANHATTAN_PLOT(grouped_regenie_ch)
    
    
  } else {
  
    regenie_file_ch2 = Channel.fromPath(params.regenie_outputs)
                             .splitText() 
                             .map { file -> tuple(file.split('/')[-1].split('\\.')[0].split('_')[1..-1].join('_'), 
                                    file)
                                  }
    region_ch = bgen_sample_ch.combine(regenie_file_ch2)
                              .filter {it[0] == it[-1].split('/')[-1].split('\\.')[0].split('_')[0].split('chr')[1]}
                                  
  }


  REGION_FILTER(region_ch)
    

  ld_input_ch = REGION_FILTER.out[0].transpose()
                .filter { it[5].name.split('_')[-2] != 'NULL' }
  				      .map { chr, pheno_id, bgen_file, bgi_file, sample_file, region_file ->
                       tuple(chr, pheno_id, region_file.name.split('_')[-2], region_file, bgen_file, bgi_file, sample_file) 
  				           }

  samples_ld_incl_ch = Channel.fromPath(params.samples_ld_incl).collect()                                  
    
    
  LD_MATRIX_CALCULATOR(ld_input_ch, samples_ld_incl_ch)
    
  
  segmentation_error_ch = LD_MATRIX_CALCULATOR.out[0].collect()
  
  
  //SEGMENTATION_FAULT_CATCHER(segmentation_error_ch)
  
    
  SUSIER(LD_MATRIX_CALCULATOR.out[0])


  coloc5_ch = SUSIER.out[0].groupTuple()
			   .combine(Channel.of('coloc5'))

  coloc3_ch = SUSIER.out[1].groupTuple()
			   .combine(Channel.of('coloc3'))

  clpp_ch = SUSIER.out[2].filter {it[1].name.split('_')[-2] != 'NULL'}
			 .groupTuple()
			 .combine(Channel.of('clpp')).view()

  
  DATA_COMBINER1(coloc5_ch)
  DATA_COMBINER2(coloc3_ch)  
  DATA_COMBINER3(clpp_ch) 
}
