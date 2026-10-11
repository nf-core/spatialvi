process IMAGE_DOWNSCALE {
    tag "${meta.id}"
    label 'process_high'

    container "community.wave.seqera.io/library/bin2cell_harmonypy_scanorama_imagecodecs_pruned:d457f4c6ebb7b325"

    input:
    tuple val(meta), path(image)
    val downscale

    output:
    tuple val(meta), path("${prefix}_downscaled.tif"), emit: image
    path "versions.yml"                              , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'downscale.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_downscaled.tif
    touch versions.yml
    """
}
