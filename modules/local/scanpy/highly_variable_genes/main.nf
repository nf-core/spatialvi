process SCANPY_HIGHLY_VARIABLE_GENES {
    tag "${meta.id}"
    label 'process_single'

    container "community.wave.seqera.io/library/harmonypy_scanorama_gcc_gxx_pruned:95f731fde0b9ddef"

    input:
    tuple val(meta), path(adata)
    val n_top_genes

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: adata
    path "versions.yml"                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}_hvg"
    flavor = task.ext.flavor ?: "seurat"
    template 'highly_variable_genes.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_hvg"
    """
    touch ${prefix}.h5ad
    touch versions.yml
    """
}
