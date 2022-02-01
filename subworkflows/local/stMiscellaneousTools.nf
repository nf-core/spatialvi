nextflow.enable.dsl=2
 
/* 
 * Include requires tasks 
 */
include { DECONVOLUTION_WITH_STDECONVOLVE   } from '../../modules/local/tasks'
include { DECONVOLUTION_WITH_SPOTLIGHT      } from '../../modules/local/tasks'
include { CLUSTERING_WITH_BAYESSPACE        } from '../../modules/local/tasks'
include { ST_SPATIALDE                      } from '../../modules/local/tasks'

/* 
 * Run miscellaneous tools
 */
workflow ST_MISCELLANEOUS_TOOLS {
 
    take:
      state
      outdir
      stRawData
      
    main:
      DECONVOLUTION_WITH_STDECONVOLVE(state, outdir, stRawData)
      DECONVOLUTION_WITH_SPOTLIGHT(state, outdir, stRawData)
      CLUSTERING_WITH_BAYESSPACE(state, outdir)
      ST_SPATIALDE(state, outdir)
       
    emit:
      ST_SPATIALDE.out
      .combine(DECONVOLUTION_WITH_STDECONVOLVE.out)
      .combine(DECONVOLUTION_WITH_SPOTLIGHT.out)
      .combine(CLUSTERING_WITH_BAYESSPACE.out)
      .collect()
     
 }
