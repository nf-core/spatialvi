//
// Subworkflow for aggregation of sample data
//

include { SDATA_MERGE                       } from "../../../modules/local/sdata/merge"
include { SCANPY_NEIGHBORS                  } from '../../../modules/local/scanpy/neighbors/main'
include { SCANPY_UMAP                       } from '../../../modules/local/scanpy/umap/main'
include { SCANPY_LEIDEN                     } from '../../../modules/local/scanpy/leiden/main'
include { SCANPY_SCANORAMA                  } from "../../../modules/local/scanpy/scanorama"
include { QUARTONOTEBOOK as INTEGRATE_SDATA } from "../../../modules/nf-core/quartonotebook/main"

workflow AGGREGATION {

    take:
    ch_sdata                       // channel: [ meta, zarr ]
    ch_adata                       // channel: [ meta, h5ad ]
    merge_sdata                    // boolean: Whether to merge sdata or not
    integrate_sdata                // boolean: Whether to integrate sdata or not
    n_neighbours                   // integer: Number of nearest neighbours to compute
    neighbours_n_pcs               // integer: Number of PCs to use for nearest neighbours
    umap_min_dist                  //   float: Minimum distance between embedded points
    umap_spread                    //   float: Scale of embedded points
    integration_cluster_resolution //   float: Integration cluster resolution
    cluster_key_added              //  string: Obs key where cluster labels are added

    main:

    // // Quarto report and extensions files
    // integration_notebook = file("${projectDir}/bin/integration.qmd", checkIfExists: true)
    // extensions = channel.fromPath("${projectDir}/assets/_extensions").collect()

    // Get sdata files only
    ch_sdata_files = ch_sdata.map { _meta, zarr -> return [zarr] }

    //
    // MODULE: Merge per-sample SpatialData objects into one
    //
    ch_merged_sdata = channel.empty()
    if (merge_sdata || integrate_sdata) {
        SDATA_MERGE (
            ch_sdata_files.collect()
        )
        ch_merged_sdata = SDATA_MERGE.out.sdata
    }

    // Conditionally run integration
    ch_integrated_sdata = channel.empty()
    ch_integrated_adata = channel.empty()
    if (integrate_sdata) {

        //
        // MODULE: Integration with Scanorama
        //
        SCANPY_SCANORAMA {
            ch_adata.map { _meta, h5ad -> h5ad }.collect()
        }
        ch_integrated = SCANPY_SCANORAMA.out.adata
            .map { h5ad -> [[id: h5ad.baseName], h5ad] }

        //
        // MODULE: Neighbourhood graph
        //
        SCANPY_NEIGHBORS (
            ch_integrated,
            n_neighbours,
            neighbours_n_pcs,
            'X_scanorama',
        )

        //
        // MODULE: UMAP
        //
        SCANPY_UMAP (
            SCANPY_NEIGHBORS.out.adata,
            umap_min_dist,
            umap_spread
        )

        //
        // MODULE: Leiden clustering
        //
        integration_cluster_key_added = cluster_key_added + "_scanorama"
        SCANPY_LEIDEN (
            SCANPY_UMAP.out.adata,
            integration_cluster_resolution,
            integration_cluster_key_added
        )
    }

    // //
    // // MODULE: Aggregate and integrate per-sample SpatialData
    // //
    // ch_integrated_sdata = channel.empty()
    // ch_integrated_adata = channel.empty()
    // if (integrate_sdata) {
    //     integration_params = [
    //         input_sdata: "merged_sdata.zarr",
    //         cluster_resolution: integration_cluster_resolution,
    //         n_hvgs: integration_n_hvgs,
    //         artifact_dir: "artifacts",
    //         output_adata: "integrated_adata.h5ad",
    //         output_sdata: "integrated_sdata.zarr"
    //     ]
    //     INTEGRATE_SDATA (
    //         [[id:"integration"], integration_notebook],
    //         integration_params,
    //         ch_merged_sdata,
    //         extensions
    //     )
    //     ch_integration_artifacts = INTEGRATE_SDATA.out.artifacts
    //         .map {
    //             _meta, artifacts ->
    //             return [artifacts]
    //         }
    //         .flatten ( )
    //         .branch { it ->
    //             adata: it[1].name.endsWith('.h5ad')
    //             sdata: it[1].name.endsWith('.zarr')
    //         }
    //     ch_integrated_adata = ch_integration_artifacts.adata
    //     ch_integrated_sdata = ch_integration_artifacts.sdata
    // }

    emit:
    merged_sdata     = ch_merged_sdata     // channel: [ zarr ]

    integrated_adata = ch_integrated_adata // channel: [ h5ad ]
    integrated_sdata = ch_integrated_sdata // channel: [ zarr ]

}
