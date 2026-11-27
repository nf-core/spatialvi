process SCANPY_LOG1P {
    tag "${meta.id}"
    label 'process_single'

    container "community.wave.seqera.io/library/bin2cell_harmonypy_scanorama_imagecodecs_pruned:d457f4c6ebb7b325"

    input:
    tuple val(meta), path(adata, stageAs: "input.h5ad", arity: '1')

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: adata
    path "versions.yml"                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'log1p.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    cat <<-END_VERSIONS > versions.yml
    "${task.process}": stub
    END_VERSIONS
    """
}
