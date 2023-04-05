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
    path(report)
    tuple val(sample_id), path(st_adata_norm, stageAs: "adata_norm.h5ad")

    output:
    tuple val(sample_id), path("*.csv"), emit: degs
    tuple val(sample_id), path("*.st_spatial_de.html")  , emit: html
    tuple val(sample_id), path("st_spatial_de_files/*")  , emit: html_files

    // path("versions.yml")               , emit: versions

    script:
    """
    quarto render "${report}" \
        --output "${sample_id}.st_spatial_de.html" \
        -P fileNameST:${st_adata_norm} \
        -P numberOfColumns:${params.st_spatial_de_ncols} \
        -P saveDEFileName:stDE.csv \
        -P saveSpatialDEFileName:st_spatial_de.csv

    mv st_spatial_de.csv "${sample_id}.st_spatial_de.csv"
    """
}
