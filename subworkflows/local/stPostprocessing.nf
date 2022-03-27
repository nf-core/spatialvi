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
      sample_ids
      outdir

    main:
      ST_CLUSTERING(  sample_ids,         outdir)
      ALL_REPORT(     ST_CLUSTERING.out,  outdir)

    emit:
      ALL_REPORT.out
 }
