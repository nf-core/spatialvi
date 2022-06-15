//
// Read ST 10x visium and SC 10x data with Scanpy and save to `anndata` file
//
process READ_ST_AND_SC_SCANPY {

    label "python_process_low"

    input:
    tuple val  (sample_id),
          path (tissue_position_list, stageAs: "SRCount/spatial/tissue_positions_list.csv"),
          path (tissue_lowres_image, stageAs: "SRCount/spatial/tissue_lowres_image.png"),
          path (tissue_hires_image, stageAs: "SRCount/spatial/tissue_hires_image.png"),
          path (scale_factors, stageAs: "SRCount/spatial/scalefactors_json.json"),
          path (barcodes, stageAs: "SRCount/raw_feature_bc_matrix/barcodes.tsv.gz"),
          path (features, stageAs: "SRCount/raw_feature_bc_matrix/features.tsv.gz"),
          path (matrix, stageAs: "SRCount/raw_feature_bc_matrix/matrix.mtx.gz")

    output:
    tuple val(sample_id), path("st_adata_raw.h5ad"), emit: st_raw
    tuple val(sample_id), path("sc_adata_raw.h5ad"), emit: sc_raw
    tuple val(sample_id), path("st_counts.npz")    , emit: st_counts
    tuple val(sample_id), path("sc_counts.npz")    , emit: sc_counts

    script:
    """
    script_read_st_data.py \
        --SRCountDir  ./SRCount \
        --outAnnData  st_adata_raw.h5ad \
        --outSTCounts st_counts.npz

    script_read_sc_data.py \
        --SRCountDir  ./SRCount \
        --outAnnData  sc_adata_raw.h5ad \
        --outSCCounts sc_counts.npz
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
// Spatial data pre-processing
//
 process ST_PREPROCESS {

    label "python_process"

    input:
    tuple val(sample_id), path(st_raw), path(st_factors)
    path(mito_data)

    output:
    tuple val(sample_id), path("*_norm.h5ad"),        emit: st_data_norm
    tuple val(sample_id), path("*_plain.h5ad"),       emit: st_data_plain
    tuple val(sample_id), path("*.st_adata_x.npz"),   emit: st_adata_x
    tuple val(sample_id), path("*.st_adata_var.npz"), emit: st_adata_var
    tuple val(sample_id), path("*.st_adata_obs.npz"), emit: st_adata_obs
    tuple val(sample_id), path("*.png"),              emit: figures

    script:
    """
    stPreprocess.py \
        --npFactorsOutputName ${st_factors} \
        --rawAdata ${st_raw} \
        --mitoFile ${mito_data} \
        --pltFigSize ${params.STpreprocess_pltFigSize} \
        --minCounts ${params.STpreprocess_minCounts} \
        --minGenes ${params.STpreprocess_minGenes} \
        --minCells ${params.STpreprocess_minCells} \
        --histplotQCmaxTotalCounts ${params.STpreprocess_histplotQCmaxTotalCounts} \
        --histplotQCminGeneCounts ${params.STpreprocess_histplotQCminGeneCounts} \
        --histplotQCbins ${params.STpreprocess_histplotQCbins} \
        --nameDataPlain ${sample_id}.st_adata_plain.h5ad \
        --nameDataNorm ${sample_id}.st_adata_norm.h5ad \
        --nameX ${sample_id}.st_adata_x.npz \
        --nameVar ${sample_id}.st_adata_var.npz \
        --nameObs ${sample_id}.st_adata_obs.npz
    """
}

//
// Single cell data pre-processing
//
 process SC_PREPROCESS {

    label "python_process"

    input:
    tuple val(sample_id), path(sc_raw), path(sc_factors)
    path(mito_data)

    output:
    tuple val(sample_id), path("*_norm.h5ad"),        emit: sc_data_norm
    tuple val(sample_id), path("*_plain.h5ad"),       emit: sc_data_plain
    tuple val(sample_id), path("*.sc_adata_x.npz"),   emit: sc_adata_x
    tuple val(sample_id), path("*.sc_adata_var.npz"), emit: sc_adata_var
    tuple val(sample_id), path("*.sc_adata_obs.npz"), emit: sc_adata_obs
    tuple val(sample_id), path("*.png"),              emit: figures

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
        --histplotQCbins=${params.SCpreprocess_histplotQCbins} \
        --nameDataPlain=${sample_id}.sc_adata_plain.h5ad \
        --nameDataNorm=${sample_id}.sc_adata_norm.h5ad \
        --nameX ${sample_id}.sc_adata_x.npz \
        --nameVar ${sample_id}.sc_adata_var.npz \
        --nameObs ${sample_id}.sc_adata_obs.npz
    """
}










//
// Spatial deconvolution with STdeconvolve
//
process DECONVOLUTION_WITH_STDECONVOLVE {

    label "r_process"

    input:
    val sample_state
    val outdir

    output:
    tuple env(sample_id), env(outpath)

    script:
    """
    characterization_STdeconvolve.R \
        --outsPath=${sample_info.st_data_dir} \
        --mtxGeneColumn=$params.STdeconvolve_mtxGeneColumn \
        --countsFactor=$params.STdeconvolve_countsFactor \
        --corpusRemoveAbove=$params.STdeconvolve_corpusRemoveAbove \
        --corpusRemoveBelow=$params.STdeconvolve_corpusRemoveBelow \
        --LDAminTopics=$params.STdeconvolve_LDAminTopics \
        --LDAmaxTopics=$params.STdeconvolve_LDAmaxTopics \
        --STdeconvolveScatterpiesSize=$params.STdeconvolve_ScatterpiesSize \
        --STdeconvolveFeaturesSizeFactor=$params.STdeconvolve_FeaturesSizeFactor
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


//
// Resolution enhancement and spatial clustering with BayesSpace
//
process CLUSTERING_WITH_BAYESSPACE {

    label "r_process"

    input:
    tuple val(sample_id), path(st_adata_x)
    tuple val(sample_id), path(st_adata_var)
    tuple val(sample_id), path(st_adata_obs)

    output:
    tuple val(sample_id), path("bayes_*.csv"), emit: tables
    tuple val(sample_id), path("*.png"), emit: figures

    script:
    """
    characterization_BayesSpace.R \
        --nameX ${st_adata_x} \
        --nameVar ${st_adata_var} \
        --nameObs ${st_adata_obs} \
        --numberHVG $params.BayesSpace_numberHVG \
        --numberPCs $params.BayesSpace_numberPCs \
        --minClusters $params.BayesSpace_minClusters \
        --maxClusters $params.BayesSpace_maxClusters \
        --optimalQ $params.BayesSpace_optimalQ \
        --STplatform $params.BayesSpace_STplatform
    """
}

//
// Clustering etc. TODO: better description
//
process ST_CLUSTERING {

    label "python_process"

    input:
    tuple val(sample_id), path(st_adata_norm), path(sc_adata_norm)

    output:
    tuple val(sample), path("*.st_*.h5ad"), emit: st_adata_processed
    tuple val(sample), path("*.sc_*.h5ad"), emit: sc_adata_processed
    path("*.png"), optional: true, emit: figures

    script:
    """
    stClusteringWorkflow.py \
        --fileNameST ${st_adata_norm} \
        --fileNameSC ${sc_adata_norm} \
        --resolution=$params.Clustering_resolution \
        --saveFileST ${sample_id}.st_adata_processed.h5ad \
        --saveFileSC ${sample_id}.sc_adata_processed.h5ad
    """
}













//
// Spatial differential expression
//
 process ST_SPATIAL_DE {

    label "python_process"

    input:
    tuple val(sample_id), path(st_data_norm)

    output:
    tuple val(sample_id), path("*.csv"), emit: degs
    path("*.png"),                       emit: figures, optional: true

    script:
    """
    stSpatialDE.py \
        --fileName=${st_data_norm} \
        --numberOfColumns=${params.SpatialDE_numberOfColumns} \
        --saveFileName=${sample_id}.stSpatialDE.csv \
        --savePlotName=${sample_id}.stSpatialDE.png
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

