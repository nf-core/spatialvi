//
// Spatial differential expression and reporting
//

include { ST_SPATIAL_DE } from '../../modules/local/st_spatial_de'
include { ST_CLUSTERING } from '../../modules/local/st_clustering'
include { REPORT_ALL    } from '../../modules/local/report_all'

workflow ST_POSTPROCESSING {

    take:
    st_data_norm

    main:
    //
    // Spatial differential expression
    //
    ST_SPATIAL_DE (
        st_data_norm
    )

    // TODO: Add reporting
    //
    // Reporting and final outputs
    //
    // REPORT_ALL     ( ST_CLUSTERING.out,  outdir)

    emit:
    spatial_degs    = ST_SPATIAL_DE.out.degs    // channel: [ val(sample), csv ]
    spatial_figures = ST_SPATIAL_DE.out.figures // channel: [ val(sample), png ]
}
