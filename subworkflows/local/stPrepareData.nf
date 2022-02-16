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
      sample_id
      outdir
      
    main:
      MITO_LOAD(sample_id, outdir)
      
    emit:
      MITO_LOAD.out

}
