//
// Calculate ST and SC sum factors for use in downstream normalisation
//
process ST_CALCULATE_SUM_FACTORS {

    tag "${sample_id}"
    label "r_process"

    input:
    tuple val(sample_id), path(st_counts)
    tuple val(sample_id), path(sc_counts)

    output:
    tuple val(sample_id), path("*.st_*.npz"), emit: st_factors
    tuple val(sample_id), path("*.sc_*.npz"), emit: sc_factors

    script:
    """
    calculateSumFactors.R \
        --npCountsOutputName=${st_counts} \
        --npFactorsOutputName=${sample_id}.st_adata_counts_in_tissue_factors.npz

    calculateSumFactors.R \
        --npCountsOutputName=${sc_counts} \
        --npFactorsOutputName=${sample_id}.sc_adata_counts_factors.npz
    """
}
