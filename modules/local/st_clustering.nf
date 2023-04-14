//
// Clustering etc.
//
process ST_CLUSTERING {

    // TODO: Add a better description
    // TODO: Add proper Conda/container directive
    // TODO: Export versions

    tag "${meta.id}"
    label "process_low"

    container "cavenel/spatialtranscriptomics"

    input:
    path(report)
    tuple val(meta), path(st_adata_norm, stageAs: "adata_norm.h5ad")

    output:
    tuple val(meta), path("st_adata_processed.h5ad"), emit: st_adata_processed
    tuple val(meta), path("st_clustering.html")     , emit: html
    tuple val(meta), path("st_clustering_files")    , emit: html_files
    // path("versions.yml")                              , emit: versions

    script:
    """
    quarto render ${report} \
        --output "st_clustering.html" \
        -P fileNameST:${st_adata_norm} \
        -P resolution:${params.st_cluster_resolution} \
        -P saveFileST:st_adata_processed.h5ad
    """
}
