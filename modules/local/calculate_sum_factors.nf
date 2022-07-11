//
// Calculate sum factors for use in downstream normalisation
//
process CALCULATE_SUM_FACTORS {

    tag "${meta.id}"
    label "process_low"

    // TODO: Add Conda/container directive

    input:
    tuple val(meta), path(counts)

    output:
    tuple val(meta), path("*.npz"), emit: factors

    script:
    """
    calculateSumFactors.R \
        --npCountsOutputName=${counts} \
        --npFactorsOutputName=${meta.id}.factors.npz
    """
}
