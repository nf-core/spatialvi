process MASK_TO_GEOJSON {
    tag "${meta.id}"
    label 'process_medium'

    container "community.wave.seqera.io/library/bin2cell_harmonypy_scanorama_imagecodecs_pruned:d457f4c6ebb7b325"

    input:
    tuple val(meta), path(mask)
    val downscale
    val object_class

    output:
    tuple val(meta), path("${prefix}_segmentation.geojson"), emit: geojson
    path "versions.yml"                                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'mask_to_geojson.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_segmentation.geojson
    touch versions.yml
    """
}
