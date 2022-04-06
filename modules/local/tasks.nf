import groovy.json.JsonSlurper

//
// Read ST 10x visium and SC 10x data with Scanpy and save to `anndata` file
//
process READ_ST_AND_SC_SCANPY {

    label "python_process_low"

    input:
    val(sample_id)
    val(outdir)

    output:
    tuple val(sample_id), path("*.st_adata_raw.h5ad"), emit: st_raw
    tuple val(sample_id), path("*.sc_adata_raw.h5ad"), emit: sc_raw
    tuple val(sample_id), path("*.st_*.npz"), emit: st_counts
    tuple val(sample_id), path("*.sc_*.npz"), emit: sc_counts

    script:
    def fileName = String.format("%s/sample_%s.json", outdir, sample_id)
    sample_info = new JsonSlurper().parse(new File(fileName))
    """
    script_read_st_data.py \
        --outsPath=${sample_info.st_data_dir} \
        --saveFile=${sample_id}.st_adata_raw.h5ad \
        --npCountsOutputName=${sample_id}.st_adata_counts_in_tissue.npz \
        --countsFile=${sample_id}.raw_feature_bc_matrix.h5 \
        --minCounts=${params.STload_minCounts} \
        --minCells=${params.STload_minCells}

    script_read_sc_data.py \
        --outsPath=${sample_info.sc_data_dir} \
        --saveFile=${sample_id}.sc_adata_raw.h5ad \
        --npCountsOutputName=${sample_id}.sc_adata_counts.npz \
        --minCounts=${params.SCload_minCounts} \
        --minCells=${params.SCload_minCells} \
        --minGenes=${params.SCload_minGenes} \
    """
}


//
// Calculate ST and SC sum factors for use in downstream normalisation
//
process ST_CALCULATE_SUM_FACTORS {

    label "r_process"

    input:
    tuple val(sample_id), path(st_counts)
    tuple val(sample_id), path(sc_counts)

    output:
    tuple val(sample_id), path("*.st_*.npz"), emit: st_factors
    tuple val(sample_id), path("*.sc_*.npz"), emit: sc_factors

    script:
    """
    calculateSumFactors.R \
        --npCountsOutputName=${st_counts} \
        --npFactorsOutputName=${sample_id}.st_adata_counts_in_tissue_factors.npz

    calculateSumFactors.R \
        --npCountsOutputName=${sc_counts} \
        --npFactorsOutputName=${sample_id}.sc_adata_counts_factors.npz
    """
}

//
// ST data preprocessing
//
 process ST_PREPROCESS {

    label "python_process"

    input:
    tuple val(sample_id), path(st_raw), path(st_factors)
    path(mito_data)

    output:
    tuple val(sample_id), path("*_plain.h5ad"), path("*_norm.h5ad"), emit: st_data
    path("*.png"), emit: figures

    script:
    """
    stPreprocess.py \
        --npFactorsOutputName=${st_factors} \
        --rawAdata=${st_raw} \
        --mitoFile=${mito_data} \
        --pltFigSize=${params.STpreprocess_pltFigSize} \
        --minCounts=${params.STpreprocess_minCounts} \
        --minGenes=${params.STpreprocess_minGenes} \
        --minCells=${params.STpreprocess_minCells} \
        --histplotQCmaxTotalCounts=${params.STpreprocess_histplotQCmaxTotalCounts} \
        --histplotQCminGeneCounts=${params.STpreprocess_histplotQCminGeneCounts} \
        --histplotQCbins=${params.STpreprocess_histplotQCbins} \
        --nameDataPlain=${sample_id}.st_adata_plain.h5ad \
        --nameDataNorm=${sample_id}.st_adata_norm.h5ad
    """
}

//
// SC data preprocessing
//
 process SC_PREPROCESS {

    label "python_process"

    input:
    tuple val(sample_id), path(sc_raw), path(sc_factors)
    path(mito_data)

    output:
    tuple val(sample_id), path("*_plain.h5ad"), path("*_norm.h5ad"), emit: sc_data
    path("*.png"), emit: figures

    script:

    """
    scPreprocess.py \
        --npFactorsOutputName=${sc_factors} \
        --rawAdata=${sc_raw} \
        --mitoFile=${mito_data} \
        --pltFigSize=${params.SCpreprocess_pltFigSize} \
        --minCounts=${params.SCpreprocess_minCounts} \
        --minGenes=${params.SCpreprocess_minGenes} \
        --minCells=${params.SCpreprocess_minCells} \
        --histplotQCmaxTotalCounts=${params.SCpreprocess_histplotQCmaxTotalCounts} \
        --histplotQCminGeneCounts=${params.SCpreprocess_histplotQCminGeneCounts} \
        --histplotQCbins=${params.SCpreprocess_histplotQCbins}
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

    Rscript $projectDir/bin/characterization_STdeconvolve.R --filePath=\${dname}/ --outsPath=${sample_info.st_data_dir} --mtxGeneColumn=$params.STdeconvolve_mtxGeneColumn --countsFactor=$params.STdeconvolve_countsFactor --corpusRemoveAbove=$params.STdeconvolve_corpusRemoveAbove --corpusRemoveBelow=$params.STdeconvolve_corpusRemoveBelow --LDAminTopics=$params.STdeconvolve_LDAminTopics --LDAmaxTopics=$params.STdeconvolve_LDAmaxTopics --STdeconvolveScatterpiesSize=$params.STdeconvolve_ScatterpiesSize --STdeconvolveFeaturesSizeFactor=$params.STdeconvolve_FeaturesSizeFactor

    if [[ -s \${dname}/STdeconvolve_prop_norm.csv ]] && \
      [[ -s \${dname}/STdeconvolve_beta_norm.csv ]] && \
      [[ -s \${dname}/STdeconvolve_sc_cluster_ids.csv ]] && \
      [[ -s \${dname}/STdeconvolve_sc_pca.csv ]] && \
      [[ -s \${dname}/STdeconvolve_sc_pca_feature_loadings.csv ]] && \
      [[ -s \${dname}/STdeconvolve_sc_cluster_markers.csv ]]
    then
      echo "completed" > "output.out" && outpath=`pwd`/output.out
    else
      echo ERROR: Output files missing. >&2
      exit 2
    fi
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

    Rscript $projectDir/bin/characterization_SPOTlight.R --filePath=\${dname}/ --outsPath=${sample_info.st_data_dir} --mtxGeneColumn=$params.SPOTlight_mtxGeneColumn --countsFactor=$params.SPOTlight_countsFactor --clusterResolution=$params.SPOTlight_clusterResolution --numberHVG=$params.SPOTlight_numberHVG --numberCellsPerCelltype=$params.SPOTlight_numberCellsPerCelltype --SPOTlightScatterpiesSize=$params.SPOTlight_ScatterpiesSize

    if [[ -s \${dname}/SPOTlight_prop_norm.csv ]] && \
      [[ -s \${dname}/SPOTlight_beta_norm.csv ]] && \
      [[ -s \${dname}/SPOTlight_sc_cluster_ids.csv ]] && \
      [[ -s \${dname}/SPOTlight_sc_pca.csv ]] && \
      [[ -s \${dname}/SPOTlight_sc_pca_feature_loadings.csv ]] && \
      [[ -s \${dname}/SPOTlight_sc_cluster_markers.csv ]]
    then
      echo "completed" > "output.out" && outpath=`pwd`/output.out
    else
      echo ERROR: Output files missing. >&2
      exit 2
    fi
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

    Rscript $projectDir/bin/characterization_BayesSpace.R --filePath=\${dname}/ --numberHVG=$params.BayesSpace_numberHVG --numberPCs=$params.BayesSpace_numberPCs --minClusters=$params.BayesSpace_minClusters --maxClusters=$params.BayesSpace_maxClusters --optimalQ=$params.BayesSpace_optimalQ --STplatform=$params.BayesSpace_STplatform

    if [[ -s \${dname}/bayes_spot_cluster.csv ]] && \
      [[ -s \${dname}/bayes_subspot_cluster_and_coord.csv ]] && \
      [[ -s \${dname}/bayes_enhanced_markers.csv ]]
    then
      echo "completed" > "output.out" && outpath=`pwd`/output.out
    else
      echo ERROR: Output files missing. >&2
      exit 2
    fi
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

    python $projectDir/bin/stSpatialDE.py --filePath=\${dname}/ --numberOfColumns=$params.SpatialDE_numberOfColumns

    if [[ -s \${dname}/stSpatialDE.csv ]]
    then
      echo "completed" > "output.out" && outpath=`pwd`/output.out
    else
      echo ERROR: Output files missing. >&2
      exit 2
    fi
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

    python $projectDir/bin/stClusteringWorkflow.py --filePath=\${dname}/ --resolution=$params.Clustering_resolution

    if [[ -s \${dname}/st_adata_processed.h5ad ]] && \
      [[ -s \${dname}/sc_adata_processed.h5ad ]]
    then
      echo "completed" > "output.out" && outpath=`pwd`/output.out
    else
      echo ERROR: Output files missing. >&2
      exit 2
    fi
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

