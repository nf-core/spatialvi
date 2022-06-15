//
// Pre-processing of spatial and single-cell data
//

include { READ_ST_AND_SC_DATA      } from '../../modules/local/read_st_and_sc_data'
include { ST_CALCULATE_SUM_FACTORS } from '../../modules/local/st_calculate_sum_factors'
include { ST_PREPROCESS            } from '../../modules/local/st_preprocess'
include { SC_PREPROCESS            } from '../../modules/local/sc_preprocess'

workflow ST_LOAD_PREPROCESS_DATA {

    take:
    ch_spatial_data

    main:

    ch_versions = Channel.empty()

    // TODO: Add file manifest or other non-hard-coded path
    //
    // Channel for mitochondrial data
    //
    ch_mito_data = Channel
        .fromPath("ftp://ftp.broadinstitute.org/distribution/metabolic/papers/Pagliarini/MitoCarta2.0/Human.MitoCarta2.0.txt")

    //
    // Read ST and SC data and save as `anndata`
    //
    READ_ST_AND_SC_DATA (
        ch_spatial_data
    )
    ch_versions = ch_versions.mix(READ_ST_AND_SC_DATA.out.versions)

    // TODO: Incorporate this step into previous one or skip it?
    //
    // Calculate sum factors used for normalisation in pre-processing
    //
    ST_CALCULATE_SUM_FACTORS (
        READ_ST_AND_SC_DATA.out.st_counts,
        READ_ST_AND_SC_DATA.out.sc_counts
    )

    //
    // Spatial pre-processing
    //
    ch_st_raw_and_factors = READ_ST_AND_SC_DATA.out.st_raw
        .join( ST_CALCULATE_SUM_FACTORS.out.st_factors )
    ST_PREPROCESS (
        ch_st_raw_and_factors,
        ch_mito_data
    )

    //
    // Single cell pre-processing
    //
    ch_sc_raw_and_factors = READ_ST_AND_SC_DATA.out.sc_raw
        .join( ST_CALCULATE_SUM_FACTORS.out.sc_factors )
    SC_PREPROCESS (
        ch_sc_raw_and_factors,
        ch_mito_data
    )

    emit:
    st_data_norm  = ST_PREPROCESS.out.st_data_norm  // channel: [ val(sample), h5ad ]
    st_data_plain = ST_PREPROCESS.out.st_data_plain // channel: [ val(sample), h5ad ]
    st_adata_x    = ST_PREPROCESS.out.st_adata_x    // channel: [ val(sample), npz ]
    st_adata_var  = ST_PREPROCESS.out.st_adata_var  // channel: [ val(sample), npz ]
    st_adata_obs  = ST_PREPROCESS.out.st_adata_obs  // channel: [ val(sample), npz ]

    sc_data_norm  = SC_PREPROCESS.out.sc_data_norm  // channel: [ val(sample), h5ad ]
    sc_data_plain = SC_PREPROCESS.out.sc_data_plain // channel: [ val(sample), h5ad ]
    sc_adata_x    = SC_PREPROCESS.out.sc_adata_x    // channel: [ val(sample), npz ]
    sc_adata_var  = SC_PREPROCESS.out.sc_adata_var  // channel: [ val(sample), npz ]
    sc_adata_obs  = SC_PREPROCESS.out.sc_adata_obs  // channel: [ val(sample), npz ]

    versions      = ch_versions                     // channel: [ version.yml ]
}
