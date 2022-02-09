nextflow.enable.dsl=2

/* 
 * Include requires tasks 
 */
include { ST_CLUSTERING   } from '../../modules/local/tasks'
include { ALL_REPORT      } from '../../modules/local/tasks'

/* 
 * Run postprocessing tools
 */
workflow ST_POSTPROCESSING {
 
    take:
      state
      sample_params
      outdir
      
    main:
      ST_CLUSTERING(state, sample_params, outdir)
      ALL_REPORT(ST_CLUSTERING.out, sample_params, outdir)
      
    emit:
      ALL_REPORT.out
 }
