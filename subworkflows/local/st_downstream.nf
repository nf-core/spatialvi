//
// Subworkflow for downstream analyses of ST data
//

include { ST_QUALITY_CONTROLS } from '../../modules/local/st_quality_controls'
include { ST_SPATIAL_DE       } from '../../modules/local/st_spatial_de'
include { ST_CLUSTERING       } from '../../modules/local/st_clustering'

workflow ST_DOWNSTREAM {

    take:
    st_sdata_raw

    main:

    ch_versions = Channel.empty()

    //
    // Report files
    //
    report_quality_controls = file("${projectDir}/bin/st_quality_controls.qmd")
    report_clustering = file("${projectDir}/bin/st_clustering.qmd")
    report_spatial_de = file("${projectDir}/bin/st_spatial_de.qmd")
    report_template = Channel.fromPath("${projectDir}/assets/_extensions").collect()

    //
    // Quality controls and filtering
    //
    ST_QUALITY_CONTROLS (
        report_quality_controls,
        report_template,
        st_sdata_raw
    )
    ch_versions = ch_versions.mix(ST_QUALITY_CONTROLS.out.versions)

    //
    // Normalisation, dimensionality reduction and clustering
    //
    ST_CLUSTERING (
        report_clustering,
        report_template,
        ST_QUALITY_CONTROLS.out.st_sdata_filtered
    )
    ch_versions = ch_versions.mix(ST_CLUSTERING.out.versions)

    //
    // Spatial differential expression
    //
    ST_SPATIAL_DE (
        report_spatial_de,
        report_template,
        ST_CLUSTERING.out.st_sdata_processed
    )
    ch_versions = ch_versions.mix(ST_SPATIAL_DE.out.versions)

    emit:
    html               = ST_QUALITY_CONTROLS.out.html               // channel: [ html ]

    st_sdata_processed = ST_CLUSTERING.out.st_sdata_processed       // channel: [ meta, h5ad]
    html               = ST_CLUSTERING.out.html                     // channel: [ html ]

    degs               = ST_SPATIAL_DE.out.degs                     // channel: [ meta, csv ]
    html               = ST_SPATIAL_DE.out.html                     // channel: [ html ]

    versions           = ch_versions                                // channel: [ versions.yml ]
}
