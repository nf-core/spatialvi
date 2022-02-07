/* 
 * Download mitochondrial genes list file
 */
 process MITO_LOAD {
    
    label "python_process"
    
    input:
    val mitoUrl
    val dataPath
     
    output:
    val("state")
              
    script:  
    """  
    #!/bin/bash

    fname=${dataPath}/`basename "${mitoUrl}"`
    echo saving to: \$fname
    
    [ ! -d ${dataPath} ] && mkdir ${dataPath}
    
    if [ ! -f \$fname ]
    then
        wget --quiet ${mitoUrl} --output-document=\$fname
    fi
    """
}


/* 
 * Read ST 10x visium and SC 10x data with scanpy and save to anndata file
 */
 process READ_ST_AND_SC_SCANPY {
 
    label "python_process"
    
    input:
    val state
    val outdir
    val stRawData
    val scRawData
    
    output:
    val outdir
         
    script:  
    """
    #!/bin/bash
      
    #sleep 5
    python $projectDir/bin/script_read_st_data.py ${stRawData} ${outdir}/st_adata_raw.h5ad raw_feature_bc_matrix.h5
    python $projectDir/bin/script_read_sc_data.py ${scRawData} ${outdir}/sc_adata_raw.h5ad
    """
}


/* 
 * Calculate ST and SC sum factors
 */
 process ST_CALCULATE_SUM_FACTORS {
    
    label "r_process"
     
    input:
    val state
    val outdir
    
    output:
    val outdir
    
    script:  
    """
    #!/bin/bash
    
    #sleep 5
    Rscript $projectDir/bin/calculateSumFactors.R ${outdir}/ st_adata_counts_in_tissue
    Rscript $projectDir/bin/calculateSumFactors.R ${outdir}/ sc_adata_counts
    """
}


/* 
 * ST data preprocessing
 */
 process ST_PREPROCESS {
    
    label "python_process"
     
    input:
    val state_factors_calculated
    val dataPath
    val mitoUrl
    val outdir
    
    output:
    val outdir
    
    script:  
    """
    #!/bin/bash
    
    mitoFile=${dataPath}/`basename "${mitoUrl}"`
     
    #sleep 5
    python $projectDir/bin/stPreprocess.py ${outdir}/ st_adata_counts_in_tissue st_adata_raw.h5ad \$mitoFile
    """
}


/* 
 * SC data preprocessing
 */
 process SC_PREPROCESS {
    
    label "python_process"
     
    input:
    val state_factors_calculated
    val dataPath
    val mitoUrl
    val outdir
    
    output:
    val outdir
    
    script:  
    """
    #!/bin/bash
    
    mitoFile=${dataPath}/`basename "${mitoUrl}"`
    
    #sleep 5
    python $projectDir/bin/scPreprocess.py ${outdir}/ sc_adata_counts sc_adata_raw.h5ad \$mitoFile
    """
}


/* 
 * ST data deconvolution with STdeconvolve
 */
 process DECONVOLUTION_WITH_STDECONVOLVE {
    
    label "r_process"
     
    input:
    val state
    val outdir
    val stRawData
    
    output:
    val outdir
    
    script:  
    """
    #!/bin/bash
      
    #sleep 5   
    Rscript $projectDir/bin/characterization_STdeconvolve.R ${outdir}/ ${stRawData}
    """
}


/* 
 * ST data deconvolution with SPOTlight
 */
 process DECONVOLUTION_WITH_SPOTLIGHT {
     
    label "r_process"
    
    input:
    val state
    val outdir
    val stRawData
    
    output:
    val outdir
    
    script:  
    """
    #!/bin/bash
    
    #sleep 5    
    Rscript $projectDir/bin/characterization_SPOTlight.R ${outdir}/ ${stRawData}
    """
}


/* 
 * Resolution enhancement and spatial clustering with BayesSpace
 */
 process CLUSTERING_WITH_BAYESSPACE {
    
    label "r_process"
     
    input:
    val state
    val outdir
    
    output:
    val outdir
    
    script:  
    """
    #!/bin/bash
       
    #sleep 5 
    Rscript $projectDir/bin/characterization_BayesSpace.R ${outdir}/
    """
}


/* 
 * SpatialDE
 */
 process ST_SPATIALDE {
     
    label "python_process"
    
    input:
    val state
    val outdir
    
    output:
    val outdir
    
    script:  
    """
    #!/bin/bash
     
    sleep 5     
    #python $projectDir/bin/stSpatialDE.py ${outdir}/ st_adata_norm.h5ad
    """
}


/* 
 * Clustering etc.
 */
 process ST_CLUSTERING {
     
    label "python_process"
    
    input:
    val state
    val outdir
    
    output:
    val outdir
    
    script:  
    """
    #!/bin/bash
     
    #sleep 5     
    python $projectDir/bin/stClusteringWorkflow.py ${outdir}/
    """
}


/* 
 * Report
 */
 process ALL_REPORT {
     
    label "python_process"
    
    input:
    val state
    val outdir
    
    output:
    val outdir
    
    script:  
    """
    #!/bin/bash
      
    #sleep 5   
    echo 1
    """
}

