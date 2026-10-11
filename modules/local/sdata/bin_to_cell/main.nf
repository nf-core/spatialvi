process SDATA_BIN_TO_CELL {
    tag "${meta.id}"
    label 'process_high'

    container "community.wave.seqera.io/library/bin2cell_harmonypy_scanorama_imagecodecs_pruned:d457f4c6ebb7b325"

    input:
    tuple val(meta), path(sdata, stageAs: "input.zarr"), path(spaceranger_dir), path(segmentation)
    val hd_bin_size

    output:
    tuple val(meta), path("${prefix}.zarr"), emit: sdata
    path "versions.yml"                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'bin_to_cell.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}.zarr
    touch ${prefix}.zarr/.zgroup
    touch versions.yml
    """
}
