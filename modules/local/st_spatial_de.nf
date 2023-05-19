//
// Spatial differential expression
//
process ST_SPATIAL_DE {

    // TODO: Add a better description
    // TODO: Find solution for Quarto with Conda

    tag "${meta.id}"
    label "process_medium"

    container "docker.io/erikfas/spatialtranscriptomics"

    input:
    path(report)
    tuple val(meta), path(st_adata_norm, stageAs: "adata_norm.h5ad")

    output:
    tuple val(meta), path("*.csv")              , emit: degs
    tuple val(meta), path("st_spatial_de.html") , emit: html
    tuple val(meta), path("st_spatial_de_files"), emit: html_files
    path("versions.yml")                        , emit: versions

    script:
    """
    quarto render ${report} \
        --output "st_spatial_de.html" \
        -P fileNameST:${st_adata_norm} \
        -P numberOfColumns:${params.st_spatial_de_ncols} \
        -P saveDEFileName:st_gde.csv \
        -P saveSpatialDEFileName:st_spatial_de.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto -v)
        leidenalg: \$(python -c "import leidenalg; print(leidenalg.__version__)")
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
        SpatialDE: \$(python -c "import SpatialDE; print(SpatialDE.__version__)")
    END_VERSIONS
    """
}
