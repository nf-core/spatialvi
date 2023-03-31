//
// Spatial differential expression and reporting
//

include { ST_SPATIAL_DE } from '../../modules/local/st_spatial_de'
include { ST_CLUSTERING } from '../../modules/local/st_clustering'
include { REPORT_ALL    } from '../../modules/local/report_all'

workflow ST_POSTPROCESSING {

    take:
    st_adata_norm

    main:

    ch_versions = Channel.empty()

    //
    // Report files
    //
    report_template_summary_clustering = file("${projectDir}/bin/stClustering.qmd")
    report_template_summary_spatial_de = file("${projectDir}/bin/stSpatialDE.qmd")

    // Clustering
    ST_CLUSTERING (
        report_template_summary_clustering,
        st_adata_norm
    )

    //
    // Spatial differential expression
    //
    ST_SPATIAL_DE (
        report_template_summary_spatial_de,
        ST_CLUSTERING.out.st_adata_processed
    )

    // TODO: Add reporting
    //
    // Reporting and final outputs
    //
    // REPORT_ALL (  )

    emit:
    spatial_degs    = ST_SPATIAL_DE.out.degs    // channel: [ val(sample), csv ]

    versions        = ch_versions               // channel: [ versions.yml ]
}
