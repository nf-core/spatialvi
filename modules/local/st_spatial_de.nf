//
// Spatial differential expression
//
 process ST_SPATIAL_DE {

    tag "${sample_id}"
    label "process_low"

    // TODO: Add Conda/container directive

    input:
    tuple val(sample_id), path(st_data_norm)

    output:
    tuple val(sample_id), path("*.csv"), emit: degs
    path("*.png")                      , emit: figures, optional: true

    script:
    """
    stSpatialDE.py \
        --fileName=${st_data_norm} \
        --numberOfColumns=${params.SpatialDE_numberOfColumns} \
        --saveFileName=${sample_id}.stSpatialDE.csv \
        --savePlotName=${sample_id}.stSpatialDE.png
    """
}
