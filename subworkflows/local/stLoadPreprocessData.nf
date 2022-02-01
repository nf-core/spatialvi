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
      state
      outdir
      stRawData
      scRawData
      dataPath
      mitoUrl
      
    main:
      READ_ST_AND_SC_SCANPY(state, outdir, stRawData, scRawData)
      ST_CALCULATE_SUM_FACTORS(READ_ST_AND_SC_SCANPY.out, outdir)
      ST_PREPROCESS(ST_CALCULATE_SUM_FACTORS.out, dataPath, mitoUrl, outdir)
      SC_PREPROCESS(ST_CALCULATE_SUM_FACTORS.out, dataPath, mitoUrl, outdir)
      
    emit:
      ST_PREPROCESS.out
      .combine(SC_PREPROCESS.out)
      .collect()
      
 }
