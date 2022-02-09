nextflow.enable.dsl=2

/* 
 * Include requires tasks 
 */
 
include { MITO_LOAD } from '../../modules/local/tasks'

/* 
 * Download and prepare data
 */
workflow ST_PREPARE_DATA {

    take:
      sample_params
      dataPath
      
    main:
      MITO_LOAD(sample_params, dataPath)
      
    emit:
      MITO_LOAD.out

}
