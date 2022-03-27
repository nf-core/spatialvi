nextflow.enable.dsl=2

/*
 * Include requires tasks
 */
include { READ_ST_AND_SC_SCANPY     } from '../../modules/local/tasks'
include { ST_CALCULATE_SUM_FACTORS  } from '../../modules/local/tasks'
include { ST_PREPROCESS             } from '../../modules/local/tasks'
include { SC_PREPROCESS             } from '../../modules/local/tasks'

/*
 * Definition of Preprocessing Workflow
 */
workflow ST_LOAD_PREPROCESS_DATA {

    take:
      sample_ids
      outdir

    main:
      READ_ST_AND_SC_SCANPY(       sample_ids,                   outdir)
      ST_CALCULATE_SUM_FACTORS(    READ_ST_AND_SC_SCANPY.out,    outdir)
      ST_PREPROCESS(               ST_CALCULATE_SUM_FACTORS.out, outdir)
      SC_PREPROCESS(               ST_CALCULATE_SUM_FACTORS.out, outdir)

    emit:
      ST_PREPROCESS.out
      .join(SC_PREPROCESS.out)

 }
