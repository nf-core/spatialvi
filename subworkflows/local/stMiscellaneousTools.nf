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

    //
    // Deconvolution with single cell data
    //
    DECONVOLUTION_WITH_STDECONVOLVE( sample_ids, outdir)
    // DECONVOLUTION_WITH_SPOTLIGHT(    sample_ids, outdir)

    //
    // Clustering etc. TODO: better description
    //
    // CLUSTERING_WITH_BAYESSPACE(      sample_ids, outdir)
    // ST_CLUSTERING ( st_data_norm.join(sc_data_norm) )

    emit:
    DECONVOLUTION_WITH_STDECONVOLVE.out
    // .join(DECONVOLUTION_WITH_STDECONVOLVE.out)
    // .join(DECONVOLUTION_WITH_SPOTLIGHT.out)
    // .join(CLUSTERING_WITH_BAYESSPACE.out)
    // st_adata_processed = ST_CLUSTERING.out.st_adata_processed // channel: [ val(sample), adata ]
    // sc_adata_processed = ST_CLUSTERING.out.sc_adata_processed // channel: [ val(sample), adata ]
    // processed_figures  = ST_CLUSTERING.out.figures            // channel: [ val(sample), png ]

 }
