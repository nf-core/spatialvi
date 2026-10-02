process SCANPY_HIGHLY_VARIABLE_GENES {
    tag "${meta.id}"
    label 'process_single'

    container "community.wave.seqera.io/library/harmonypy_scanorama_gcc_gxx_pruned:95f731fde0b9ddef"

    input:
    tuple val(meta), path(h5ad, stageAs: "input.h5ad", arity: '1')
    val n_hvgs
    val flavor
    val write_adata

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: adata, optional: true
    tuple val(meta), path("var/*.pkl")     , emit: var
    tuple val(meta), path("uns/*.pkl")     , emit: uns
    path "versions.yml"                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'highly_variable_genes.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def touch_adata = write_adata.toString() == 'true' ? "touch ${prefix}.h5ad" : ''
    """
    ${touch_adata}
    mkdir -p var uns
    touch var/highly_variable_genes.pkl
    touch uns/hvg.pkl
    touch versions.yml
    """
}
