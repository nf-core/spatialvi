//
// Test helper: give the last gene the name of the first gene, creating a
// duplicate in `var_names`; used to test modules that require unique names.
//
// TODO: Should be replaced with a proper test dataset that contains duplicate
// gene names, once the overall structure and format of Python modules is set.
//
process DUPLICATE_VAR_NAMES {
    tag "${meta.id}"
    label 'process_single'

    container "community.wave.seqera.io/library/harmonypy_scanorama_gcc_gxx_pruned:95f731fde0b9ddef"

    input:
    tuple val(meta), path(h5ad, stageAs: "input.h5ad", arity: '1')

    output:
    tuple val(meta), path("${prefix}.h5ad"), emit: adata

    script:
    prefix = task.ext.prefix ?: "${meta.id}_duplicated"
    """
    #!/usr/bin/env python3
    import anndata as ad

    adata = ad.read_h5ad("${h5ad}")
    adata.var_names = [*adata.var_names[:-1], adata.var_names[0]]
    adata.write_h5ad("${prefix}.h5ad")
    """
}
