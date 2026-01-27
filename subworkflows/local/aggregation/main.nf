//
// Subworkflow for aggregation of sample data
//

include { MERGE_SDATA                       } from "../../../modules/local/merge_sdata"
include { QUARTONOTEBOOK as INTEGRATE_SDATA } from "../../../modules/nf-core/quartonotebook/main"

workflow AGGREGATION {

    take:
    ch_sdata                       // channel: [ meta, zarr ]
    merge_sdata                    // boolean: Whether to merge sdata or not
    integrate_sdata                // boolean: Whether to integrate sdata or not
    integration_cluster_resolution // float  : Integration cluster resolution
    integration_n_hvgs             // integer: Number of HVGs to use for integration

    main:

    // Quarto report and extensions files
    integration_notebook = file("${projectDir}/bin/integration.qmd", checkIfExists: true)
    extensions = channel.fromPath("${projectDir}/assets/_extensions").collect()

    // Get sdata files only
    ch_sdata_files = ch_sdata
        .map {
            _meta, zarr ->
            return [zarr]
        }

    //
    // MODULE: Merge per-sample SpatialData objects into one
    //
    ch_merged_sdata = channel.empty()
    if (merge_sdata || integrate_sdata) {
        MERGE_SDATA (
            ch_sdata_files.collect()
        )
        ch_merged_sdata = MERGE_SDATA.out.sdata
    }

    //
    // MODULE: Aggregate and integrate per-sample SpatialData
    //
    ch_integrated_sdata = channel.empty()
    ch_integrated_adata = channel.empty()
    if (integrate_sdata) {
        integration_params = [
            input_sdata: "merged_sdata.zarr",
            cluster_resolution: integration_cluster_resolution,
            n_hvgs: integration_n_hvgs,
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
        ch_integration_artifacts = INTEGRATE_SDATA.out.artifacts
            .map {
                _meta, artifacts ->
                return [artifacts]
            }
            .flatten ( )
            .branch { it ->
                adata: it[1].name.endsWith('.h5ad')
                sdata: it[1].name.endsWith('.zarr')
            }
        ch_integrated_adata = ch_integration_artifacts.adata
        ch_integrated_sdata = ch_integration_artifacts.sdata
    }

    emit:
    merged_sdata     = ch_merged_sdata     // channel: [ zarr ]

    integrated_adata = ch_integrated_adata // channel: [ h5ad ]
    integrated_sdata = ch_integrated_sdata // channel: [ zarr ]

}
