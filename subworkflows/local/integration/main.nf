//
// Subworkflow for aggregation of sample data
//

include { SCANPY_NEIGHBORS                                     } from '../../../modules/local/scanpy/neighbors/main'
include { SCANPY_UMAP                                          } from '../../../modules/local/scanpy/umap/main'
include { SCANPY_LEIDEN                                        } from '../../../modules/local/scanpy/leiden/main'
include { SCANPY_SCANORAMA                                     } from "../../../modules/local/scanpy/scanorama"
include { SDATA_UPDATE_TABLE as SDATA_UPDATE_TABLE_INTEGRATION } from '../../../modules/local/sdata/update_table/main'
include { QUARTONOTEBOOK as REPORT_INTEGRATED                  } from "../../../modules/nf-core/quartonotebook/main"

workflow INTEGRATION {

    take:
    ch_sdata_merged                // channel: [ meta, zarr ]
    ch_adata                       // channel: [ meta, h5ad ]
    integration_method             //  string: Integration method to use
    n_neighbours                   // integer: Number of nearest neighbours to compute
    neighbours_n_pcs               // integer: Number of PCs to use for nearest neighbours
    umap_min_dist                  //   float: Minimum distance between embedded points
    umap_spread                    //   float: Scale of embedded points
    integration_cluster_resolution //   float: Integration cluster resolution

    main:

    //
    // MODULE: Integration with Scanorama
    //
    if (integration_method == 'scanorama') {
        SCANPY_SCANORAMA {
            ch_adata.map { _meta, h5ad -> h5ad }.collect()
        }
        ch_adata_integrated = SCANPY_SCANORAMA.out.adata
    }

    // Add `meta` to integrated AnnData channel
    ch_integrated = ch_adata_integrated
        .map { h5ad -> [[id: integration_method], h5ad] }

    //
    // MODULE: Neighbourhood graph
    //
    use_rep = 'X_' + integration_method
    SCANPY_NEIGHBORS (
        ch_integrated,
        n_neighbours,
        neighbours_n_pcs,
        use_rep
    )

    //
    // MODULE: UMAP
    //
    umap_key_added = 'X_umap_' + integration_method
    SCANPY_UMAP (
        SCANPY_NEIGHBORS.out.adata,
        umap_min_dist,
        umap_spread,
        umap_key_added
    )

    //
    // MODULE: Leiden clustering
    //
    leiden_key_added = 'clusters_' + integration_method
    SCANPY_LEIDEN (
        SCANPY_UMAP.out.adata,
        integration_cluster_resolution,
        leiden_key_added
    )
    ch_adata_integrated = SCANPY_LEIDEN.out.adata
        .map { _meta, h5ad -> h5ad }

    //
    // MODULE: Update SpatialData tables
    //
    library_key = 'library_id'
    ch_sdata_adata = ch_sdata_merged
        .map  { zarr -> [[id: integration_method], zarr]}
        .join ( SCANPY_LEIDEN.out.adata )
    SDATA_UPDATE_TABLE_INTEGRATION (
        ch_sdata_adata,
        library_key
    )
    ch_sdata_integrated = SDATA_UPDATE_TABLE_INTEGRATION.out.sdata

    //
    // MODULE: Aggregate and integrate per-sample SpatialData
    //
    integration_notebook = file("${projectDir}/bin/report-integrated.qmd", checkIfExists: true)
    extensions = channel.fromPath("${projectDir}/assets/_extensions").collect()
    integration_params = [
        input_adata: integration_method + ".h5ad",
        input_sdata: integration_method + ".zarr",
        artifact_dir: "artifacts",
    ]
    integration_inputs = ch_sdata_integrated
        .map { _meta, zarr -> zarr }
        .mix ( ch_adata_integrated )
        .collect()
    REPORT_INTEGRATED (
        [[id: integration_method], integration_notebook],
        integration_params,
        integration_inputs,
        extensions
    )

    emit:
    sdata            = ch_sdata_integrated // channel: [ zarr ]
    sdata_merged     = ch_sdata_merged     // channel: [ zarr ]
    adata            = ch_adata_integrated // channel: [ h5ad ]
}
