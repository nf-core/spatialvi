//
// Subworkflow for aggregation of sample data
//

include { SDATA_MERGE                                          } from "../../../modules/local/sdata/merge"
include { SCANPY_NEIGHBORS                                     } from '../../../modules/local/scanpy/neighbors/main'
include { SCANPY_UMAP                                          } from '../../../modules/local/scanpy/umap/main'
include { SCANPY_LEIDEN                                        } from '../../../modules/local/scanpy/leiden/main'
include { SCANPY_SCANORAMA                                     } from "../../../modules/local/scanpy/scanorama"
include { SDATA_UPDATE_TABLE as SDATA_UPDATE_TABLE_INTEGRATION } from '../../../modules/local/sdata/update_table/main'
include { QUARTONOTEBOOK as REPORT_INTEGRATED                  } from "../../../modules/nf-core/quartonotebook/main"

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

    // Get sdata files only
    ch_sdata_files = ch_sdata.map { _meta, zarr -> return [zarr] }

    //
    // MODULE: Merge per-sample SpatialData objects into one
    //
    ch_sdata_merged = channel.empty()
    if (merge_sdata || integrate_sdata) {
        SDATA_MERGE (
            ch_sdata_files.collect()
        )
        ch_sdata_merged = SDATA_MERGE.out.sdata
    }

    // Conditionally run integration
    ch_sdata_integrated = channel.empty()
    ch_adata_integrated = channel.empty()
    if (integrate_sdata) {

        //
        // MODULE: Integration with Scanorama
        //
        SCANPY_SCANORAMA {
            ch_adata.map { _meta, h5ad -> h5ad }.collect()
        }
        ch_integrated = SCANPY_SCANORAMA.out.adata
            .map { h5ad -> [[id: "integrated"], h5ad] }

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
        ch_adata_integrated = SCANPY_LEIDEN.out.adata
            .map { _meta, h5ad -> h5ad }

        //
        // MODULE: Update SpatialData tables
        //
        ch_sdata_adata = ch_sdata_merged
            .map  { zarr -> [[id: "integrated"], zarr]}
            .join ( SCANPY_LEIDEN.out.adata )
        SDATA_UPDATE_TABLE_INTEGRATION (
            ch_sdata_adata,
            'library_id'
        )
        ch_sdata_integrated = SDATA_UPDATE_TABLE_INTEGRATION.out.sdata

        //
        // MODULE: Aggregate and integrate per-sample SpatialData
        //
        integration_notebook = file("${projectDir}/bin/report-integrated.qmd", checkIfExists: true)
        extensions = channel.fromPath("${projectDir}/assets/_extensions").collect()
        integration_params = [
            input_adata: "integrated.h5ad",
            input_sdata: "integrated.zarr",
            artifact_dir: "artifacts",
        ]
        integration_inputs = ch_sdata_integrated
            .map { _meta, zarr -> zarr }
            .mix ( ch_adata_integrated )
            .collect()
        REPORT_INTEGRATED (
            [[id: "integrated"], integration_notebook],
            integration_params,
            integration_inputs,
            extensions
        )

    }

    emit:
    sdata            = ch_sdata_integrated // channel: [ zarr ]
    sdata_merged     = ch_sdata_merged     // channel: [ zarr ]
    adata            = ch_adata_integrated // channel: [ h5ad ]
}
