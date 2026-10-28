process SQUIDPY_SPATIAL_NEIGHBORS {
    tag "${meta.id}"
    label 'process_medium'

    container "community.wave.seqera.io/library/harmonypy_scanorama_gcc_gxx_pruned:95f731fde0b9ddef"

    input:
    tuple val(meta), path(adata, stageAs: "input.h5ad", arity: '1')
    val coord_type
    val n_neighs
    val write_adata

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: adata, optional: true
    tuple val(meta), path("obsp/*.pkl")    , emit: obsp
    tuple val(meta), path("uns/*.pkl")     , emit: uns
    path "versions.yml"                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'spatial_neighbors.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def touch_adata = write_adata.toString() == 'true' ? "touch ${prefix}.h5ad" : ''
    """
    ${touch_adata}
    mkdir -p obsp uns
    touch obsp/spatial_connectivities.pkl
    touch obsp/spatial_distances.pkl
    touch uns/spatial_neighbors.pkl
    touch versions.yml
    """
}
