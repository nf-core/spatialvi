//
// Pre-processing of ST data
//

include { CALCULATE_SUM_FACTORS } from '../../modules/local/calculate_sum_factors'
include { ST_PREPROCESS         } from '../../modules/local/st_preprocess'

workflow PREPROCESS_ST_DATA {

    take:
    st_counts
    st_raw
    mito_data

    main:

    ch_versions = Channel.empty()

    // TODO: Incorporate this step into `READ_ST_AND_SC_DATA` or skip it?
    //
    // Calculate sum factors used for normalisation in pre-processing
    //
    CALCULATE_SUM_FACTORS (
        st_counts
    )
    // ch_versions = ch_versions.mix(CALCULATE_SUM_FACTORS.out.versions)

    //
    // Spatial pre-processing
    //
    ST_PREPROCESS (
        st_raw.join(CALCULATE_SUM_FACTORS.out.factors),
        mito_data
    )
    // ch_versions = ch_versions.mix(ST_PREPROCESS.out.versions)

    emit:
    st_data_norm  = ST_PREPROCESS.out.st_data_norm  // channel: [ val(sample), h5ad ]
    st_data_plain = ST_PREPROCESS.out.st_data_plain // channel: [ val(sample), h5ad ]
    st_adata_x    = ST_PREPROCESS.out.st_adata_x    // channel: [ val(sample), npz ]
    st_adata_var  = ST_PREPROCESS.out.st_adata_var  // channel: [ val(sample), npz ]
    st_adata_obs  = ST_PREPROCESS.out.st_adata_obs  // channel: [ val(sample), npz ]

    versions      = ch_versions                     // channel: [ version.yml ]
}
