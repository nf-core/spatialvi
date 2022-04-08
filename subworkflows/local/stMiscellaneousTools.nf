nextflow.enable.dsl=2

/*
 * Include requires tasks
 */
include { DECONVOLUTION_WITH_STDECONVOLVE } from '../../modules/local/tasks'
include { DECONVOLUTION_WITH_SPOTLIGHT    } from '../../modules/local/tasks'
include { CLUSTERING_WITH_BAYESSPACE      } from '../../modules/local/tasks'

/*
 * Run miscellaneous tools
 */
workflow ST_MISCELLANEOUS_TOOLS {

    take:
    sample_ids
    outdir

    main:
    DECONVOLUTION_WITH_STDECONVOLVE( sample_ids, outdir)
    // DECONVOLUTION_WITH_SPOTLIGHT(    sample_ids, outdir)
    // CLUSTERING_WITH_BAYESSPACE(      sample_ids, outdir)

    emit:
    DECONVOLUTION_WITH_STDECONVOLVE.out
    // .join(DECONVOLUTION_WITH_STDECONVOLVE.out)
    // .join(DECONVOLUTION_WITH_SPOTLIGHT.out)
    // .join(CLUSTERING_WITH_BAYESSPACE.out)

 }
