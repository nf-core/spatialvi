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
      sample_params
      outdir
      
    main:
      DECONVOLUTION_WITH_STDECONVOLVE(state, sample_params, outdir)
      DECONVOLUTION_WITH_SPOTLIGHT(state, sample_params, outdir)
      CLUSTERING_WITH_BAYESSPACE(state, sample_params, outdir)
      ST_SPATIALDE(state, sample_params, outdir)
       
    emit:
      ST_SPATIALDE.out
      .join(DECONVOLUTION_WITH_STDECONVOLVE.out)
      .join(DECONVOLUTION_WITH_SPOTLIGHT.out)
      .join(CLUSTERING_WITH_BAYESSPACE.out)
     
 }
