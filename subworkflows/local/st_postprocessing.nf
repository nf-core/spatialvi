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
    report_template_summary = file("${projectDir}/bin/stClustering.qmd")

    // Clustering
    ST_CLUSTERING (
        report_template_summary,
        st_adata_norm
    )

    //
    // Spatial differential expression
    //
    ST_SPATIAL_DE (
        st_adata_norm
    )

    // TODO: Add reporting
    //
    // Reporting and final outputs
    //
    // REPORT_ALL (  )

    emit:
    spatial_degs    = ST_SPATIAL_DE.out.degs    // channel: [ val(sample), csv ]
    spatial_figures = ST_SPATIAL_DE.out.figures // channel: [ val(sample), png ]

    versions        = ch_versions               // channel: [ versions.yml ]
}
