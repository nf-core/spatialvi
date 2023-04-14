//
// Spatial differential expression
//
process ST_SPATIAL_DE {

    // TODO: Add proper Conda/container directive
    // TODO: Export versions

    tag "${meta.id}"
    label "process_medium"

    container "cavenel/spatialtranscriptomics"

    input:
    path(report)
    tuple val(meta), path(st_adata_norm, stageAs: "adata_norm.h5ad")

    output:
    tuple val(meta), path("*/*.csv"), emit: degs
    tuple val(meta), path("*/st_spatial_de.html")  , emit: html

    // path("versions.yml")               , emit: versions

    script:
    """
    quarto render ${report} \
        --output "st_spatial_de.html" \
        -P fileNameST:${st_adata_norm} \
        -P numberOfColumns:${params.st_spatial_de_ncols} \
        -P saveDEFileName:st_gde.csv \
        -P saveSpatialDEFileName:st_spatial_de.csv

    mkdir -p ${meta.id}
    # mv st_gde.csv ${meta.id}/st_gde.csv
    mv st_spatial_de.csv ${meta.id}/st_spatial_de.csv
    mv st_spatial_de.html ${meta.id}/st_spatial_de.html
    """
}
