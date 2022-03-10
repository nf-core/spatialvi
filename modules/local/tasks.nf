import groovy.json.JsonSlurper

/* 
 * Download mitochondrial genes list file
 */
 process MITO_LOAD {
    
    label "python_process_low"
    
    input:
    val sample_id
    val outdir
     
    output:
    tuple val(sample_id), env(outpath)
              
    script:
    def fileName = String.format("%s/sample_%s.json", outdir, sample_id)
    sample_info = new JsonSlurper().parse(new File(fileName))
        
    """  
    #!/bin/bash
    
    mitoUrl="ftp://ftp.broadinstitute.org/distribution/metabolic/papers/Pagliarini/MitoCarta2.0/${sample_info.species}.MitoCarta2.0.txt"

    fname=${outdir}/`basename "\${mitoUrl}"`
    echo saving to: \$fname
    
    [ ! -d ${outdir} ] && mkdir ${outdir}
    
    if [ ! -f \$fname ]
    then
        wget --quiet \${mitoUrl} --output-document=\$fname
    fi

    # Output is the mito file
    echo "completed" > "output.out" && outpath=`pwd`/output.out
    """
}


/* 
 * Read ST 10x visium and SC 10x data with scanpy and save to anndata file
 */
 process READ_ST_AND_SC_SCANPY {
 
    label "python_process_low"
    
    input:
    tuple val(sample_id), val(state)
    val outdir
    
    output:
    tuple val(sample_id), env(outpath)
         
    script:
    def fileName = String.format("%s/sample_%s.json", outdir, sample_id)
    sample_info = new JsonSlurper().parse(new File(fileName))
      
    """
    #!/bin/bash
    
    dname=${outdir}/${sample_id}
      
    [ ! -d \${dname} ] && mkdir \${dname}

    python $projectDir/bin/script_read_st_data.py --outsPath=${sample_info.st_data_dir} --saveFile=\${dname}/st_adata_raw.h5ad --countsFile=raw_feature_bc_matrix.h5 --npCountsOutputName=st_adata_counts_in_tissue.npz

    python $projectDir/bin/script_read_sc_data.py --outsPath=${sample_info.sc_data_dir} --saveFile=\${dname}/sc_adata_raw.h5ad --npCountsOutputName=sc_adata_counts.npz

    # Outputs from script 1: saveFile, npCountsOutputName
    # Outputs from script 2: saveFile, npCountsOutputName
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
    tuple val(sample_id), val(state)
    val outdir
    
    output:
    tuple val(sample_id), env(outpath)
    
    script:  
    """
    #!/bin/bash
    
    dname=${outdir}/${sample_id}
    
    Rscript $projectDir/bin/calculateSumFactors.R --filePath=\${dname}/ --npCountsOutputName=st_adata_counts_in_tissue.npz --npFactorsOutputName=st_adata_counts_in_tissue_factors.npz
    Rscript $projectDir/bin/calculateSumFactors.R --filePath=\${dname}/ --npCountsOutputName=sc_adata_counts.npz --npFactorsOutputName=sc_adata_counts_factors.npz

    # Outputs from script 1: factors file
    # Outputs from script 2: factors file
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
    tuple val(sample_id), val(state)
    val outdir
    
    output:
    tuple val(sample_id), env(outpath)
    
    script:
    def fileName = String.format("%s/sample_%s.json", outdir, sample_id)
    sample_info = new JsonSlurper().parse(new File(fileName))
    
    """
    #!/bin/bash
    
    dname=${outdir}/${sample_id}
    
    mitoFile=${outdir}/${sample_info.species}.MitoCarta2.0.txt
    
    python $projectDir/bin/stPreprocess.py --filePath=\${dname}/ --npFactorsOutputName=st_adata_counts_in_tissue.npz --rawAdata=st_adata_raw.h5ad --mitoFile=\$mitoFile

    # Outputs from script: plots & normalized h5ad, X, Var, Obs
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
    tuple val(sample_id), val(state)
    val outdir
    
    output:
    tuple val(sample_id), env(outpath)
    
    script:
    def fileName = String.format("%s/sample_%s.json", outdir, sample_id)
    sample_info = new JsonSlurper().parse(new File(fileName))
      
    """
    #!/bin/bash
    
    dname=${outdir}/${sample_id}
    
    mitoFile=${outdir}/${sample_info.species}.MitoCarta2.0.txt
    
    python $projectDir/bin/scPreprocess.py --filePath=\${dname}/ --npFactorsOutputName=sc_adata_counts.npz --rawAdatasc_adata_raw.h5ad --mitoFile=\$mitoFile

    # Outputs from script: plots & normalized h5ad, X, Var, Obs
    echo "completed" > "output.out" && outpath=`pwd`/output.out
    """
}












/* 
 * ST data deconvolution with STdeconvolve
 */
 process DECONVOLUTION_WITH_STDECONVOLVE {
    
    label "r_process"
     
    input:
    val sample_state
    val outdir
    
    output:
    tuple env(sample_id), env(outpath)
    
    script:
    def sample_id_gr = sample_state[0]
    def fileName = String.format("%s/sample_%s.json", outdir, sample_id_gr)
    sample_info = new JsonSlurper().parse(new File(fileName))
        
    """
    #!/bin/bash
    
    sample_id=${sample_id_gr}
    
    dname=${outdir}/\${sample_id}
       
    Rscript $projectDir/bin/characterization_STdeconvolve.R --filePath=\${dname}/ --outsPath=${sample_info.st_data_dir}

    # Output from script: plots & 6 files: prop, beta, and sc files
    echo "completed" > "output.out" && outpath=`pwd`/output.out
    """
}


/* 
 * ST data deconvolution with SPOTlight
 */
 process DECONVOLUTION_WITH_SPOTLIGHT {
     
    label "r_process"
    
    input:
    val sample_state
    val outdir
    
    output:
    tuple env(sample_id), env(outpath)
    
    script:
    def sample_id_gr = sample_state[0]
    def fileName = String.format("%s/sample_%s.json", outdir, sample_id_gr)
    sample_info = new JsonSlurper().parse(new File(fileName))
        
    """
    #!/bin/bash
    
    sample_id=${sample_id_gr}
    
    dname=${outdir}/\${sample_id}
        
    Rscript $projectDir/bin/characterization_SPOTlight.R --filePath=\${dname}/ --outsPath=${sample_info.st_data_dir}

    # Output from script: plots & 6 files: prop, beta, and sc files
    echo "completed" > "output.out" && outpath=`pwd`/output.out
    """
}


/* 
 * Resolution enhancement and spatial clustering with BayesSpace
 */
 process CLUSTERING_WITH_BAYESSPACE {
    
    label "r_process"
     
    input:
    val sample_state
    val outdir
    
    output:
    tuple env(sample_id), env(outpath)
    
    script:
    def sample_id_gr = sample_state[0]
    def fileName = String.format("%s/sample_%s.json", outdir, sample_id_gr)
    sample_info = new JsonSlurper().parse(new File(fileName))
        
    """
    #!/bin/bash
    
    sample_id=${sample_id_gr}
    
    dname=${outdir}/\${sample_id}
    
    Rscript $projectDir/bin/characterization_BayesSpace.R --filePath=\${dname}/

    # Output from script: plots & 3 files: clusters, subclusters, enhanced genes
    echo "completed" > "output.out" && outpath=`pwd`/output.out
    """
}


/* 
 * SpatialDE
 */
 process ST_SPATIALDE {
     
    label "python_process"
    
    input:
    val sample_state
    val outdir
    
    output:
    tuple env(sample_id), env(outpath)
    
    script:
    def sample_id_gr = sample_state[0]
    def fileName = String.format("%s/sample_%s.json", outdir, sample_id_gr)
    sample_info = new JsonSlurper().parse(new File(fileName))
        
    """
    #!/bin/bash
    
    sample_id=${sample_id_gr}
    
    dname=${outdir}/\${sample_id}
       
    python $projectDir/bin/stSpatialDE.py --filePath=\${dname}/

    # Output from script: plots & ...
    echo "completed" > "output.out" && outpath=`pwd`/output.out
    """
}
















/* 
 * Clustering etc.
 */
 process ST_CLUSTERING {
     
    label "python_process"
    
    input:
    val sample_state
    val outdir
    
    output:
    tuple env(sample_id), env(outpath)
    
    script:
    def sample_id_gr = sample_state[0]
    def fileName = String.format("%s/sample_%s.json", outdir, sample_id_gr)
    sample_info = new JsonSlurper().parse(new File(fileName))
        
    """
    #!/bin/bash
    
    sample_id=${sample_id_gr}
    
    dname=${outdir}/\${sample_id}
         
    python $projectDir/bin/stClusteringWorkflow.py --filePath=\${dname}/

    # Output from script: ...
    echo "completed" > "output.out" && outpath=`pwd`/output.out
    """
}


/* 
 * Report
 */
 process ALL_REPORT {
     
    label "python_process_low"
    
    input:
    val sample_state
    val outdir
    
    output:
    tuple env(sample_id), env(outpath)
    
    script:
    def sample_id_gr = sample_state[0]
    def fileName = String.format("%s/sample_%s.json", outdir, sample_id_gr)
    sample_info = new JsonSlurper().parse(new File(fileName))
        
    """
    #!/bin/bash
    
    sample_id=${sample_id_gr}
    
    dname=${outdir}/\${sample_id}
    
    echo \${dname}/  
    echo "completed" > "output.out" && outpath=`pwd`/output.out
    """
}

