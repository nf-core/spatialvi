//
// Pre-processing of ST data
//

include { ST_PREPROCESS         } from '../../modules/local/st_preprocess'

workflow PREPROCESS_ST_DATA {

    take:
    st_raw
    mito_data

    main:

    ch_versions = Channel.empty()

    //
    // Report files
    //
    report_template_summary = file("${projectDir}/bin/stPreprocess.qmd")

    //
    // Spatial pre-processing
    //
    ST_PREPROCESS (
        report_template_summary,
        st_raw,
        mito_data
    )
    // ch_versions = ch_versions.mix(ST_PREPROCESS.out.versions)

    emit:
    st_data_norm  = ST_PREPROCESS.out.st_data_norm  // channel: [ val(sample), h5ad ]
    st_data_plain = ST_PREPROCESS.out.st_data_plain // channel: [ val(sample), h5ad ]

    versions      = ch_versions                     // channel: [ version.yml ]
}
