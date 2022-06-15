//
// Read ST 10x visium and SC 10x data with Scanpy and save to `anndata` file
//
process READ_ST_AND_SC_DATA {

    tag "${sample_id}"
    label "process_low"

    conda (params.enable_conda ? "conda-forge::scanpy=1.7.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scanpy:1.7.2--pyhdfd78af_0' :
        'quay.io/biocontainers/scanpy:1.7.2--pyhdfd78af_0' }"

    input:
    tuple val  (sample_id),
          path (tissue_position_list, stageAs: "SRCount/spatial/tissue_positions_list.csv"),
          path (tissue_lowres_image , stageAs: "SRCount/spatial/tissue_lowres_image.png"),
          path (tissue_hires_image  , stageAs: "SRCount/spatial/tissue_hires_image.png"),
          path (scale_factors       , stageAs: "SRCount/spatial/scalefactors_json.json"),
          path (barcodes            , stageAs: "SRCount/raw_feature_bc_matrix/barcodes.tsv.gz"),
          path (features            , stageAs: "SRCount/raw_feature_bc_matrix/features.tsv.gz"),
          path (matrix              , stageAs: "SRCount/raw_feature_bc_matrix/matrix.mtx.gz")

    output:
    tuple val(sample_id), path("st_adata_raw.h5ad"), emit: st_raw
    tuple val(sample_id), path("sc_adata_raw.h5ad"), emit: sc_raw
    tuple val(sample_id), path("st_counts.npz")    , emit: st_counts
    tuple val(sample_id), path("sc_counts.npz")    , emit: sc_counts
    path("versions.yml")                           , emit: versions

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

     cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """
}
