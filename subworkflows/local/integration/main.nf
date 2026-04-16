//
// Subworkflow for aggregation of sample data
//

include { ADATA_MERGE                                          } from "../../../modules/local/adata/merge"
include { QUARTONOTEBOOK as REPORT_INTEGRATED                  } from "../../../modules/nf-core/quartonotebook"
include { SCANPY_HARMONY                                       } from "../../../modules/local/scanpy/harmony"
include { SCANPY_SCANORAMA                                     } from "../../../modules/local/scanpy/scanorama"
include { SDATA_UPDATE_TABLE as SDATA_UPDATE_TABLE_INTEGRATION } from "../../../modules/local/sdata/update_table"

include { CLUSTERING                                           } from "../../../subworkflows/local/clustering"

workflow INTEGRATION {

    take:
    ch_sdata_merged                // channel: [ meta, zarr ]
    ch_adata                       // channel: [ meta, h5ad ]
    integration_method             //  string: Integration method to use
    n_neighbors                   // integer: Number of nearest neighbors to compute
    neighbors_n_pcs               // integer: Number of PCs to use for nearest neighbors
    umap_min_dist                  //   float: Minimum distance between embedded points
    umap_spread                    //   float: Scale of embedded points
    integration_cluster_resolution //   float: Integration cluster resolution

    main:

    //
    // MODULE: Merge AnnData objects
    //
    ch_adata_collected = ch_adata
        .toSortedList { a, b -> a[0].id <=> b[0].id }
        .flatMap()
        .map { _meta, zarr -> zarr }
        .collect()
    ADATA_MERGE (
        ch_adata_collected,
        'inner',      // join
        'library_id', // label
        'true',       // preserve_var
        'true'        // preserve_spatial
    )
    ch_adata_merged = ADATA_MERGE.out.adata
        .map { h5ad -> [[id: integration_method], h5ad] }

    //
    // MODULE: Integration
    //
    if (integration_method == 'harmony') {
        SCANPY_HARMONY (
            ch_adata_merged,
            'library_id', // key
            'X_harmony'   // adjusted_basis
        )
        ch_adata_integrated = SCANPY_HARMONY.out.adata
    } else if (integration_method == 'scanorama') {
        SCANPY_SCANORAMA (
            ch_adata_merged,
            'library_id', // key
            'X_scanorama' // embedding_added
        )
        ch_adata_integrated = SCANPY_SCANORAMA.out.adata
    }

    //
    // SUBWORKFLOW: Clustering
    //
    use_rep = 'X_' + integration_method
    umap_key_added = 'X_umap_' + integration_method
    leiden_key_added = 'clusters_' + integration_method
    CLUSTERING (
        ch_adata_integrated,
        n_neighbors,
        neighbors_n_pcs,
        use_rep,
        umap_min_dist,
        umap_spread,
        umap_key_added,
        integration_cluster_resolution,
        leiden_key_added
    )
    ch_adata_clustered = CLUSTERING.out.adata
        .map { _meta, h5ad -> h5ad }

    //
    // MODULE: Update SpatialData tables
    //
    library_key = 'library_id'
    ch_sdata_adata = ch_sdata_merged
        .map  { zarr -> [[id: integration_method], zarr]}
        .join ( CLUSTERING.out.adata )
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
        sample_col: 'library_id',
        cluster_col: 'clusters_' + integration_method,
        artifact_dir: "artifacts",
    ]
    integration_inputs = ch_sdata_integrated
        .map { _meta, zarr -> zarr }
        .mix ( ch_adata_clustered )
        .collect()
    REPORT_INTEGRATED (
        [[id: integration_method], integration_notebook],
        integration_params,
        integration_inputs,
        extensions
    )

    emit:
    sdata       = ch_sdata_integrated               // channel: [ zarr ]
    adata       = ch_adata_clustered                // channel: [ h5ad ]

    html        = REPORT_INTEGRATED.out.html        // channel: [ meta, html ]
    notebook    = REPORT_INTEGRATED.out.notebook    // channel: [ meta, qmd ]
    params_yaml = REPORT_INTEGRATED.out.params_yaml // channel: [ meta, yml ]
    artifacts   = REPORT_INTEGRATED.out.artifacts   // channel: [ meta, dir ]
}
