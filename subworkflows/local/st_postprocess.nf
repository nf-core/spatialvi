//
// Spatial differential expression and reporting
//

include { ST_SPATIAL_DE } from '../../modules/local/st_spatial_de'
include { ST_CLUSTERING } from '../../modules/local/st_clustering'

workflow ST_POSTPROCESS {

    take:
    st_adata_norm

    main:

    ch_versions = Channel.empty()

    //
    // Report files
    //
    report_clustering = file("${projectDir}/bin/st_clustering.qmd")
    report_spatial_de = file("${projectDir}/bin/st_spatial_de.qmd")

    //
    // Clustering
    //
    ST_CLUSTERING (
        report_clustering,
        st_adata_norm
    )
    ch_versions = ch_versions.mix(ST_CLUSTERING.out.versions)

    //
    // Spatial differential expression
    //
    ST_SPATIAL_DE (
        report_spatial_de,
        ST_CLUSTERING.out.st_adata_processed
    )
    ch_versions = ch_versions.mix(ST_SPATIAL_DE.out.versions)

    emit:
    spatial_degs    = ST_SPATIAL_DE.out.degs    // channel: [ val(sample), csv ]

    versions        = ch_versions               // channel: [ versions.yml ]
}
