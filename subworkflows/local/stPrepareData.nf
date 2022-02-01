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
      species
      mitoUrl
      dataPath
      
    main:
      MITO_LOAD(mitoUrl, dataPath)
      
    emit:
      MITO_LOAD.out

}
