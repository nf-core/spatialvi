//
// Spatial differential expression
//
process ST_SPATIAL_DE {

    // TODO: Add proper Conda/container directive
    // TODO: Export versions

    tag "${sample_id}"
    label "process_low"

    container "cavenel/spatialtranscriptomics"

    input:
    path(report_template_summary)
    tuple val(sample_id), path(st_adata_norm, stageAs: "adata_norm.h5ad")

    output:
    tuple val(sample_id), path("*.csv"), emit: degs
    tuple val(sample_id), path("*.stSpatialDE.html")  , emit: report
    tuple val(sample_id), path("stSpatialDE_files/*")  , emit: report_files

    // path("versions.yml")               , emit: versions

    script:
    """
    quarto render "${report_template_summary}" --output "${sample_id}.stSpatialDE.html" \
        -P fileNameST:${st_adata_norm} \
        -P numberOfColumns:${params.SpatialDE_numberOfColumns} \
        -P saveDEFileName:stDE.csv \
        -P saveSpatialDEFileName:stSpatialDE.csv

    # mv stDE.csv "${sample_id}.stDE.csv"
    mv stSpatialDE.csv "${sample_id}.stSpatialDE.csv"
    """
}
