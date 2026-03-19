process SDATA_MERGE {
    label 'process_low'

    container "community.wave.seqera.io/library/harmonypy_scanorama_gcc_gxx_pruned:95f731fde0b9ddef"

    input:
    path(sdata, stageAs: "?/*")

    output:
    path("sdata_merged.zarr"), emit: sdata
    path "versions.yml"      , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "sdata_merged"
    template 'merge.py'

    stub:
    prefix = task.ext.prefix ?: "sdata_merged"
    """
    touch ${prefix}.zarr
    touch versions.yml
    """
}
