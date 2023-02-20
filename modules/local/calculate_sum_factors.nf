//
// Calculate sum factors for use in downstream normalisation
//
process CALCULATE_SUM_FACTORS {

    // TODO: Add proper Conda/container directive
    // Export versions

    tag "${meta.id}"
    label "process_low"

    container "erikfas/spatialtranscriptomics"

    input:
    tuple val(meta), path(counts)

    output:
    tuple val(meta), path("*.npz"), emit: factors
    // path("versions.yml")          , emit: versions

    script:
    """
    calculateSumFactors.R \
        --npCountsOutputName=${counts} \
        --npFactorsOutputName=${meta.id}.factors.npz
    """
}
