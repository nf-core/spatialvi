//
// Spatial differential expression and reporting
//

include { ST_SPATIAL_DE } from '../../modules/local/st_spatial_de'
include { ST_CLUSTERING } from '../../modules/local/st_clustering'

workflow ST_POSTPROCESS {

    take:
    st_adata_filtered

    main:

    ch_versions = Channel.empty()

    //
    // Report files
    //
    report_clustering = file("${projectDir}/bin/st_clustering.qmd")
    report_spatial_de = file("${projectDir}/bin/st_spatial_de.qmd")
    report_template = Channel.fromPath("${projectDir}/assets/_extensions")

    //
    // Normalisation, dimensionality reduction and clustering
    //
    ST_CLUSTERING (
        report_clustering,
        report_template,
        st_adata_filtered
    )
    ch_versions = ch_versions.mix(ST_CLUSTERING.out.versions)

    //
    // Spatial differential expression
    //
    ST_SPATIAL_DE (
        report_spatial_de,
        report_template,
        ST_CLUSTERING.out.st_adata_processed
    )
    ch_versions = ch_versions.mix(ST_SPATIAL_DE.out.versions)

    emit:
    st_adata_processed = ST_CLUSTERING.out.st_adata_processed // channel: [ meta, h5ad]
    html               = ST_CLUSTERING.out.html               // channel: [ html ]

    degs               = ST_SPATIAL_DE.out.degs               // channel: [ meta, csv ]
    html               = ST_SPATIAL_DE.out.html               // channel: [ html ]

    versions           = ch_versions                          // channel: [ versions.yml ]
}
