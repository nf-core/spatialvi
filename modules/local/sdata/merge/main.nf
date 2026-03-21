process SDATA_MERGE {
    label 'process_low'

    container "community.wave.seqera.io/library/harmonypy_scanorama_gcc_gxx_pruned:95f731fde0b9ddef"

    input:
    path(sdata, stageAs: "?/*")

    output:
    path("${prefix}.zarr"), emit: sdata
    path "versions.yml"   , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "merged"
    template 'merge.py'

    stub:
    prefix = task.ext.prefix ?: "merged"
    """
    touch ${prefix}.zarr
    touch versions.yml
    """
}
