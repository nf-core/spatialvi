//
// Pre-processing of ST data
//

include { ST_QUALITY_CONTROLS } from '../../modules/local/st_quality_controls'

workflow ST_PREPROCESS {

    take:
    st_raw

    main:

    ch_versions = Channel.empty()

    //
    // Report files
    //
    report = file("${projectDir}/bin/st_quality_controls.qmd")
    report_template = Channel.fromPath("${projectDir}/assets/_extensions")

    //
    // Spatial pre-processing
    //
    ST_QUALITY_CONTROLS (
        report,
        report_template,
        st_raw
    )
    ch_versions = ch_versions.mix(ST_QUALITY_CONTROLS.out.versions)

    emit:
    st_data_norm  = ST_QUALITY_CONTROLS.out.st_data_norm  // channel: [ val(sample), h5ad ]

    versions      = ch_versions                               // channel: [ version.yml ]
}
