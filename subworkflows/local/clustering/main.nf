include { SCANPY_LEIDEN                } from '../../../modules/local/scanpy/leiden/main'
include { SCANPY_NEIGHBORS             } from '../../../modules/local/scanpy/neighbors/main'
include { SCANPY_UMAP                  } from '../../../modules/local/scanpy/umap/main'

workflow CLUSTERING {

    take:
    ch_adata           // channel: [ meta, h5ad ]
    n_neighbours       // integer: Number of nearest neighbours to compute
    neighbours_n_pcs   // integer: Number of PCs to use for nearest neighbours
    use_rep            //  string: Embedding to use for nearest neighbours
    umap_min_dist      //   float: Minimum distance between embedded points
    umap_spread        //   float: Scale of embedded points
    umap_key_added     //  string: Key in .obsm for storing UMAP coordinates
    cluster_resolution //   float: Spot clustering resolution
    leiden_key_added   //  string: Key in .obs for storing cluster labels

    main:

    //
    // MODULE: Neighbourhood graph
    //
    SCANPY_NEIGHBORS (
        ch_adata,
        n_neighbours,
        neighbours_n_pcs,
        use_rep
    )

    //
    // MODULE: UMAP
    //
    SCANPY_UMAP (
        SCANPY_NEIGHBORS.out.adata,
        umap_min_dist,
        umap_spread,
        umap_key_added
    )

    //
    // MODULE: Leiden clustering
    //
    SCANPY_LEIDEN (
        SCANPY_UMAP.out.adata,
        cluster_resolution,
        leiden_key_added
    )

    emit:
    adata = SCANPY_LEIDEN.out.adata // channel: [ meta, h5ad ]
}
