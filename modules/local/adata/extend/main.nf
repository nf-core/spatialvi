process ADATA_EXTEND {
    tag "${meta.id}"
    label 'process_single'

    container "community.wave.seqera.io/library/harmonypy_scanorama_gcc_gxx_pruned:95f731fde0b9ddef"

    input:
    tuple (
        val(meta),
        path(h5ad,   stageAs: "base.h5ad", arity: '1'),
        path(obs,    stageAs: "obs/*"),
        path(var,    stageAs: "var/*"),
        path(obsm,   stageAs: "obsm/*"),
        path(varm,   stageAs: "varm/*"),
        path(obsp,   stageAs: "obsp/*"),
        path(uns,    stageAs: "uns/*"),
        path(layers, stageAs: "layers/*")
    )
    val align
    val overwrite

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: adata
    path "versions.yml"                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'extend.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.h5ad
    touch versions.yml
    """
}
