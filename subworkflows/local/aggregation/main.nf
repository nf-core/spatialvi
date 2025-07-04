//
// Subworkflow for aggregation of sample data
//

include { MERGE_SDATA                       } from "../../../modules/local/merge_sdata"
include { QUARTONOTEBOOK as INTEGRATE_SDATA } from "../../../modules/nf-core/quartonotebook/main"

workflow AGGREGATION {

    take:
    ch_sdata // Channel: [ meta, zarr ]

    main:

    ch_versions = Channel.empty()

    // Quarto report and extensions files
    integration_notebook = file("${projectDir}/bin/integration.qmd", checkIfExists: true)
    extensions = Channel.fromPath("${projectDir}/assets/_extensions").collect()

    // Get sdata files only
    ch_sdata_files = ch_sdata
        | map {
            meta, zarr ->
            return [zarr]
        }

    //
    // MODULE: Merge per-sample SpatialData objects into one
    //
    ch_merged_sdata = Channel.empty()
    if (params.merge_sdata || params.integrate_sdata) {
        MERGE_SDATA (
            ch_sdata_files.collect()
        )
        ch_versions = ch_versions.mix(MERGE_SDATA.out.versions)
        ch_merged_sdata = MERGE_SDATA.out.sdata
    }

    //
    // MODULE: Aggregate and integrate per-sample SpatialData
    //
    ch_integrated_sdata = Channel.empty()
    ch_integrated_adata = Channel.empty()
    if (params.integrate_sdata) {
        integration_params = [
            input_sdata: "merged_sdata.zarr",
            cluster_resolution: params.integration_cluster_resolution,
            n_hvgs: params.integration_n_hvgs,
            artifact_dir: "artifacts",
            output_adata: "integrated_adata.h5ad",
            output_sdata: "integrated_sdata.zarr"
        ]
        INTEGRATE_SDATA (
            [[id:"integration"], integration_notebook],
            integration_params,
            ch_merged_sdata,
            extensions
        )
        ch_versions = ch_versions.mix(INTEGRATE_SDATA.out.versions)
        ch_integration_artifacts = INTEGRATE_SDATA.out.artifacts
            | map {
                meta, artifacts ->
                return [artifacts]
            }
            | flatten()
            | branch {
                adata: it[1].name.endsWith('.h5ad')
                sdata: it[1].name.endsWith('.zarr')
            }
        ch_integrated_adata = ch_integration_artifacts.adata
        ch_integrated_sdata = ch_integration_artifacts.sdata
    }

    emit:
    merged_sdata     = ch_merged_sdata     // channel: [ aggregated-sdata.zarr ]

    integrated_adata = ch_integrated_adata // channel: [ integrated_adata.h5ad ]
    integrated_sdata = ch_integrated_sdata // channel: [ integrated_sdata.zarr ]

    versions         = ch_versions         // channel: [ versions.yml ]

}
