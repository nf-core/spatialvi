process SDATA_UPDATE_TABLE {
    tag "${meta.id}"
    label 'process_single'

    container "community.wave.seqera.io/library/harmonypy_scanorama_gcc_gxx_pruned:95f731fde0b9ddef"

    input:
    tuple val(meta), path(sdata), path(adata)

    output:
    tuple val(meta), path("${prefix}_updated.zarr"), emit: sdata
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'update_table.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_updated.zarr
    touch ${prefix}_updated.zarr/.zgroup
    touch versions.yml
    """
}