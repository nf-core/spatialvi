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

    //
    // Read ST and SC data and save as `anndata`
    //
    READ_ST_AND_SC_SCANPY ( sample_ids, outdir )

    //
    // Calculate sum factors used for normalisation in pre-processing
    //
    ST_CALCULATE_SUM_FACTORS( READ_ST_AND_SC_SCANPY.out.st_counts,
                              READ_ST_AND_SC_SCANPY.out.sc_counts)

    ST_PREPROCESS( ST_CALCULATE_SUM_FACTORS.out.st_factors, outdir)
    SC_PREPROCESS( ST_CALCULATE_SUM_FACTORS.out.sc_factors, outdir)

    emit:
    ST_PREPROCESS.out.join(SC_PREPROCESS.out)
}
