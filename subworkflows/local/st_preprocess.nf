//
// Pre-processing of ST data
//

include { ST_QC_AND_NORMALISATION } from '../../modules/local/st_qc_and_normalisation'

workflow ST_PREPROCESS {

    take:
    st_raw

    main:

    ch_versions = Channel.empty()

    //
    // Report files
    //
    report = file("${projectDir}/bin/st_qc_and_normalisation.qmd")

    //
    // Spatial pre-processing
    //
    ST_QC_AND_NORMALISATION (
        report,
        st_raw
    )
    ch_versions = ch_versions.mix(ST_QC_AND_NORMALISATION.out.versions)

    emit:
    st_data_norm  = ST_QC_AND_NORMALISATION.out.st_data_norm  // channel: [ val(sample), h5ad ]
    st_data_plain = ST_QC_AND_NORMALISATION.out.st_data_plain // channel: [ val(sample), h5ad ]

    versions      = ch_versions                               // channel: [ version.yml ]
}
