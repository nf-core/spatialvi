include { SQUIDPY_INTERACTION_MATRIX } from "../../../modules/local/squidpy/interaction_matrix"
include { SQUIDPY_NHOOD_ENRICHMENT   } from "../../../modules/local/squidpy/nhood_enrichment"
include { SQUIDPY_SPATIAL_AUTOCORR   } from "../../../modules/local/squidpy/spatial_autocorr"
include { SQUIDPY_SPATIAL_NEIGHBORS  } from "../../../modules/local/squidpy/spatial_neighbors"

workflow SPATIAL {

    take:
    ch_adata            // channel: [ meta, h5ad ]
    spatial_coord_type  //  string: Type of spatial coordinate system
    spatial_n_neighbors // integer: Number of spatial neighbors to use
    svg_autocorr_method //  string: Spatial autocorrelation method to use

    main:

    //
    // Spatial neighbors
    //
    SQUIDPY_SPATIAL_NEIGHBORS (
        ch_adata,
        spatial_coord_type,
        spatial_n_neighbors
    )

    //
    // Spatial neighborhood enrichment analysis
    //
    cluster_key = 'clusters'
    SQUIDPY_NHOOD_ENRICHMENT (
        SQUIDPY_SPATIAL_NEIGHBORS.out.adata,
        cluster_key
    )

    //
    // Spatial interaction matrix
    //
    SQUIDPY_INTERACTION_MATRIX (
        SQUIDPY_NHOOD_ENRICHMENT.out.adata,
        cluster_key
    )

    //
    // Spatial autocorrelation (spatially variable genes)
    //
    SQUIDPY_SPATIAL_AUTOCORR (
        SQUIDPY_INTERACTION_MATRIX.out.adata,
        svg_autocorr_method
    )

    emit:
    adata   = SQUIDPY_SPATIAL_AUTOCORR.out.adata // channel: [ meta, h5ad ]
    svg_csv = SQUIDPY_SPATIAL_AUTOCORR.out.csv   // channel: [ meta, csv ]
}
