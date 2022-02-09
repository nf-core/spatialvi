/* 
 * Download mitochondrial genes list file
 */
 process MITO_LOAD {
    
    label "python_process_low"
    
    input:
    tuple val(sample_id), val(species), val(st_data_dir), val(sc_data_dir)
    val dataPath
     
    output:
    tuple val(sample_id), env(outpath)
              
    script:  
    """  
    #!/bin/bash
    
    mitoUrl="ftp://ftp.broadinstitute.org/distribution/metabolic/papers/Pagliarini/MitoCarta2.0/${species}.MitoCarta2.0.txt"

    fname=${dataPath}/`basename "\${mitoUrl}"`
    echo saving to: \$fname
    
    [ ! -d ${dataPath} ] && mkdir ${dataPath}
    
    if [ ! -f \$fname ]
    then
        wget --quiet \${mitoUrl} --output-document=\$fname
    fi
    echo "completed" > "output.out" && outpath=`pwd`/output.out
    """
}


/* 
 * Read ST 10x visium and SC 10x data with scanpy and save to anndata file
 */
 process READ_ST_AND_SC_SCANPY {
 
    label "python_process_low"
    
    input:
    tuple val(s), path(p)
    tuple val(sample_id), val(species), val(st_data_dir), val(sc_data_dir)
    val outdir
    
    output:
    tuple val(sample_id), env(outpath)
         
    script:  
    """
    #!/bin/bash
    
    dname=${outdir}/${sample_id}
      
    [ ! -d \${dname} ] && mkdir \${dname}
    
    python $projectDir/bin/script_read_st_data.py ${st_data_dir} \${dname}/st_adata_raw.h5ad raw_feature_bc_matrix.h5
    python $projectDir/bin/script_read_sc_data.py ${sc_data_dir} \${dname}/sc_adata_raw.h5ad
    echo "completed" > "output.out" && outpath=`pwd`/output.out
    """
}


/* 
 * Calculate ST and SC sum factors
 */
 process ST_CALCULATE_SUM_FACTORS {
    
    label "r_process"
    cpus 2
    memory '8.GB'
     
    input:
    tuple val(s), path(p)
    tuple val(sample_id), val(species), val(st_data_dir), val(sc_data_dir)
    val outdir
    
    output:
    tuple val(sample_id), env(outpath)
    
    script:  
    """
    #!/bin/bash
    
    dname=${outdir}/${sample_id}
    
    Rscript $projectDir/bin/calculateSumFactors.R \${dname}/ st_adata_counts_in_tissue
    Rscript $projectDir/bin/calculateSumFactors.R \${dname}/ sc_adata_counts
    echo "completed" > "output.out" && outpath=`pwd`/output.out
    """
}


/* 
 * ST data preprocessing
 */
 process ST_PREPROCESS {
    
    label "python_process"
    cpus 2
    memory '4.GB'
     
    input:
    tuple val(s), path(p)
    tuple val(sample_id), val(species), val(st_data_dir), val(sc_data_dir)
    val dataPath
    val outdir
    
    output:
    tuple val(sample_id), env(outpath)
    
    script:  
    """
    #!/bin/bash
    
    dname=${outdir}/${sample_id}
    
    mitoFile=${dataPath}/${species}.MitoCarta2.0.txt
    
    python $projectDir/bin/stPreprocess.py \${dname}/ st_adata_counts_in_tissue st_adata_raw.h5ad \$mitoFile
    echo "completed" > "output.out" && outpath=`pwd`/output.out
    """
}


/* 
 * SC data preprocessing
 */
 process SC_PREPROCESS {
    
    label "python_process"
    cpus 2
    memory '4.GB'
         
    input:
    tuple val(s), path(p)
    tuple val(sample_id), val(species), val(st_data_dir), val(sc_data_dir)
    val dataPath
    val outdir
    
    output:
    tuple val(sample_id), env(outpath)
    
    script:  
    """
    #!/bin/bash
    
    dname=${outdir}/${sample_id}
    
    mitoFile=${dataPath}/${species}.MitoCarta2.0.txt
    
    python $projectDir/bin/scPreprocess.py \${dname}/ sc_adata_counts sc_adata_raw.h5ad \$mitoFile
    echo "completed" > "output.out" && outpath=`pwd`/output.out
    """
}


/* 
 * ST data deconvolution with STdeconvolve
 */
 process DECONVOLUTION_WITH_STDECONVOLVE {
    
    label "r_process"
     
    input:
    tuple val(s), path(p)
    tuple val(sample_id), val(species), val(st_data_dir), val(sc_data_dir)
    val outdir
    
    output:
    tuple val(sample_id), env(outpath)
    
    script:  
    """
    #!/bin/bash
    
    dname=${outdir}/${sample_id}
       
    Rscript $projectDir/bin/characterization_STdeconvolve.R \${dname}/ ${st_data_dir}
    echo "completed" > "output.out" && outpath=`pwd`/output.out
    """
}


/* 
 * ST data deconvolution with SPOTlight
 */
 process DECONVOLUTION_WITH_SPOTLIGHT {
     
    label "r_process"
    
    input:
    tuple val(s), path(p)
    tuple val(sample_id), val(species), val(st_data_dir), val(sc_data_dir)
    val outdir
    
    output:
    tuple val(sample_id), env(outpath)
    
    script:  
    """
    #!/bin/bash
    
    dname=${outdir}/${sample_id}
        
    Rscript $projectDir/bin/characterization_SPOTlight.R \${dname}/ ${st_data_dir}
    echo "completed" > "output.out" && outpath=`pwd`/output.out
    """
}


/* 
 * Resolution enhancement and spatial clustering with BayesSpace
 */
 process CLUSTERING_WITH_BAYESSPACE {
    
    label "r_process"
     
    input:
    tuple val(s), path(p)
    tuple val(sample_id), val(species), val(st_data_dir), val(sc_data_dir)
    val outdir
    
    output:
    tuple val(sample_id), env(outpath)
    
    script:  
    """
    #!/bin/bash
    
    dname=${outdir}/${sample_id}
    
    Rscript $projectDir/bin/characterization_BayesSpace.R \${dname}/
    echo "completed" > "output.out" && outpath=`pwd`/output.out
    """
}


/* 
 * SpatialDE
 */
 process ST_SPATIALDE {
     
    label "python_process"
    
    input:
    tuple val(s), path(p)
    tuple val(sample_id), val(species), val(st_data_dir), val(sc_data_dir)
    val outdir
    
    output:
    tuple val(sample_id), env(outpath)
    
    script:  
    """
    #!/bin/bash
    
    dname=${outdir}/${sample_id}
       
    python $projectDir/bin/stSpatialDE.py \${dname}/ st_adata_norm.h5ad
    echo "completed" > "output.out" && outpath=`pwd`/output.out
    """
}


/* 
 * Clustering etc.
 */
 process ST_CLUSTERING {
     
    label "python_process"
    
    input:
    tuple val(s), path(p)
    tuple val(sample_id), val(species), val(st_data_dir), val(sc_data_dir)
    val outdir
    
    output:
    tuple val(sample_id), env(outpath)
    
    script:  
    """
    #!/bin/bash
    
    dname=${outdir}/${sample_id}
         
    python $projectDir/bin/stClusteringWorkflow.py \${dname}/
    echo "completed" > "output.out" && outpath=`pwd`/output.out
    """
}


/* 
 * Report
 */
 process ALL_REPORT {
     
    label "python_process_low"
    
    input:
    tuple val(s), path(p)
    tuple val(sample_id), val(species), val(st_data_dir), val(sc_data_dir)
    val outdir
    
    output:
    tuple val(sample_id), env(outpath)
    
    script:  
    """
    #!/bin/bash
    
    dname=${outdir}/${sample_id}
    
    echo \${dname}/  
    echo "completed" > "output.out" && outpath=`pwd`/output.out
    """
}

