//
// Read ST 10x visium and SC 10x data with Scanpy and save to `anndata` file
//
process ST_READ_DATA {

    tag "${meta.id}"
    label "process_low"

    conda "conda-forge::scanpy=1.7.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scanpy:1.7.2--pyhdfd78af_0' :
        'quay.io/biocontainers/scanpy:1.7.2--pyhdfd78af_0' }"

    input:
    tuple val (meta),
        path (tissue_positions_list, stageAs: "SRCount/spatial/tissue_positions_list.csv"),
        path (tissue_lowres_image  , stageAs: "SRCount/spatial/tissue_lowres_image.png"),
        path (tissue_hires_image   , stageAs: "SRCount/spatial/tissue_hires_image.png"),
        path (scale_factors        , stageAs: "SRCount/spatial/scalefactors_json.json"),
        path (barcodes             , stageAs: "SRCount/raw_feature_bc_matrix/barcodes.tsv.gz"),
        path (features             , stageAs: "SRCount/raw_feature_bc_matrix/features.tsv.gz"),
        path (matrix               , stageAs: "SRCount/raw_feature_bc_matrix/matrix.mtx.gz")

    output:
    tuple val(meta), path("st_adata_raw.h5ad"), emit: st_raw
    path("versions.yml")                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    read_st_data.py \
        --SRCountDir ./SRCount \
        --outAnnData st_adata_raw.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """
}
