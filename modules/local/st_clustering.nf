//
// Clustering etc.
//
process ST_CLUSTERING {

    // TODO: Add a better description
    // TODO: Add proper Conda/container directive
    // TODO: Export versions

    label "process_medium"

    container "cavenel/spatialtranscriptomics"

    input:
    path(report_template_summary)
    tuple val(sample_id), path(st_adata_norm, stageAs: "adata_norm.h5ad")

    output:
    tuple val(sample_id), path("*.st_*.h5ad"), emit: st_adata_processed
    tuple val(sample_id), path("*.stClustering.html")  , emit: report
    tuple val(sample_id), path("stClustering_files/*")  , emit: report_files
    // path("versions.yml")                  , emit: versions

    script:
    """
    quarto render "${report_template_summary}" --output "${sample_id}.stClustering.html" \
        -P fileNameST:${st_adata_norm} \
        -P resolution:$params.Clustering_resolution \
        -P saveFileST:st_adata_processed.h5ad

    mv st_adata_processed.h5ad "${sample_id}.st_adata_processed.h5ad"
    """
}
