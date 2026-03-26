process SDATA_UPDATE_TABLE {
    tag "${meta.id}"
    label 'process_single'

    container "community.wave.seqera.io/library/harmonypy_scanorama_gcc_gxx_pruned:95f731fde0b9ddef"

    input:
    tuple val(meta), path(sdata, stageAs: "input.zarr"), path(adata)
    val library_key

    output:
    tuple val(meta), path("${prefix}.zarr"), emit: sdata
    tuple val(meta), path(adata)           , emit: adata
    path "versions.yml"                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'update_table.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}.zarr
    touch ${prefix}.zarr/.zgroup
    touch versions.yml
    """
}
