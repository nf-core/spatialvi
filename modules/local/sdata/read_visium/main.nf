process SDATA_READ_VISIUM {
    tag "${meta.id}"
    label 'process_low'

    container "community.wave.seqera.io/library/harmonypy_scanorama_gcc_gxx_pruned:95f731fde0b9ddef"

    input:
    tuple val(meta), path(spaceranger_dir, arity: '1')
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
