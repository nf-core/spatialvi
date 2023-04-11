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
    tuple val(sample_id), path("*/*.csv"), emit: degs
    tuple val(sample_id), path("*/st_spatial_de.html")  , emit: html

    // path("versions.yml")               , emit: versions

    script:
    """
    quarto render "${report}" \
        --output "st_spatial_de.html" \
        -P fileNameST:${st_adata_norm} \
        -P numberOfColumns:${params.st_spatial_de_ncols} \
        -P saveDEFileName:st_gde.csv \
        -P saveSpatialDEFileName:st_spatial_de.csv

    mkdir "${sample_id}" -p
    
    # mv st_gde.csv "${sample_id}/st_gde.csv"
    mv st_spatial_de.csv "${sample_id}/st_spatial_de.csv"
    mv st_spatial_de.html "${sample_id}/st_spatial_de.html"
    """
}
