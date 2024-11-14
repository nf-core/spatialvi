//
// Subworkflow for aggregation of sample data
//

include { MERGE_SDATA } from '../../modules/local/merge_sdata'

workflow AGGREGATION {

    take:
    ch_sdata // Channel: [ meta, zarr ]

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Merge per-sample SpatialData objects into one
    //
    ch_sdata_files = ch_sdata
        | map {
            meta, zarr ->
            return [zarr]
        }
    MERGE_SDATA (
        ch_sdata_files.collect()
    )
    ch_versions = ch_versions.mix(MERGE_SDATA.out.versions)
    ch_merged_sdata = MERGE_SDATA.out.sdata

    emit:
    merged_sdata = ch_merged_sdata // channel: [ aggregated-sdata.zarr ]
    versions     = ch_versions     // channel: [ versions.yml ]

}
