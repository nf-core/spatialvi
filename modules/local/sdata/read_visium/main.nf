process SDATA_READ_VISIUM {
    tag "${meta.id}"
    label 'process_low'

    container "community.wave.seqera.io/library/bin2cell_harmonypy_scanorama_imagecodecs_pruned:d457f4c6ebb7b325"

    input:
    tuple val(meta), path(spaceranger_dir, arity: '1'), val(read_cells)
    val(hd_bin_size)

    output:
    tuple val(meta), path("${prefix}.zarr"), emit: sdata
    path "versions.yml"                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "The SDATA_READ_VISIUM module does not support Conda/Mamba, please use Docker / Singularity / Podman instead."
    }
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'read_visium.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}.zarr
    touch ${prefix}.zarr/.zgroup
    touch versions.yml
    """
}
