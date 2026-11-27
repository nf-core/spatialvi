process SDATA_MERGE {
    label 'process_low'

    container "community.wave.seqera.io/library/bin2cell_harmonypy_scanorama_imagecodecs_pruned:d457f4c6ebb7b325"

    input:
    path(sdata, stageAs: "?/*", arity: '1..*')

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
    mkdir -p ${prefix}.zarr
    touch ${prefix}.zarr/.zgroup
    cat <<-END_VERSIONS > versions.yml
    "${task.process}": stub
    END_VERSIONS
    """
}
