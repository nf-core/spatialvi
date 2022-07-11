//
// Pre-processing of SC data
//

include { CALCULATE_SUM_FACTORS } from '../../modules/local/calculate_sum_factors'
include { SC_PREPROCESS         } from '../../modules/local/sc_preprocess'

workflow PREPROCESS_SC_DATA {

    take:
    sc_counts
    sc_raw
    mito_data

    main:

    ch_versions = Channel.empty()

    // TODO: Incorporate this step into `READ_ST_AND_SC_DATA` or skip it?
    //
    // Calculate sum factors used for normalisation in pre-processing
    //
    CALCULATE_SUM_FACTORS (
        sc_counts
    )
    // ch_versions = ch_versions.mix(CALCULATE_SUM_FACTORS.out.versions)

    //
    // Spatial pre-processing
    //
    SC_PREPROCESS (
        sc_raw.join(CALCULATE_SUM_FACTORS.out.factors),
        mito_data
    )
    // ch_versions = ch_versions.mix(SC_PREPROCESS.out.versions)

    emit:
    sc_data_norm  = SC_PREPROCESS.out.sc_data_norm  // channel: [ val(sample), h5ad ]
    sc_data_plain = SC_PREPROCESS.out.sc_data_plain // channel: [ val(sample), h5ad ]
    sc_adata_x    = SC_PREPROCESS.out.sc_adata_x    // channel: [ val(sample), npz ]
    sc_adata_var  = SC_PREPROCESS.out.sc_adata_var  // channel: [ val(sample), npz ]
    sc_adata_obs  = SC_PREPROCESS.out.sc_adata_obs  // channel: [ val(sample), npz ]

    versions      = ch_versions                     // channel: [ version.yml ]
}
