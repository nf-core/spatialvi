#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/spatialvi
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/spatialvi
    Website: https://nf-co.re/spatialvi
    Slack  : https://nfcore.slack.com/channels/spatialvi
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SPATIALVI               } from "./workflows/spatialvi"
include { PIPELINE_INITIALISATION } from "./subworkflows/local/utils_nfcore_spatialvi_pipeline"
include { PIPELINE_COMPLETION     } from "./subworkflows/local/utils_nfcore_spatialvi_pipeline"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_SPATIALVI {

    take:
    samplesheet                    //    file: samplesheet read in from --input
    spaceranger_reference          //     dir: /path/to/reference
    spaceranger_probeset           //    file: /path/to/csv
    hd_bin_size                    // integer: Bin size for Visium HD
    skip_downstream                // boolean: Whether to skip downstream steps or not
    qc_min_counts                  // integer: Minimum UMIs per spot
    qc_min_genes                   // integer: Minimum genes per spot
    qc_min_spots                   // integer: Minimum spots per gene
    qc_mito_threshold              //   float: Maximum mito. content per spot
    qc_ribo_threshold              //   float: Minimum ribo. content per spot
    qc_hb_threshold                //   float: Maximum haem. content per spot
    normalize_target_sum           //  string: Target sum of total count normalization
    n_highly_variable_genes        // integer: Number of HVGs to use
    hvg_flavor                     //  string: Flavor for HVG calculations
    n_principal_components         // integer: Number of principal components to compute
    pca_use_highly_variable        // boolean: Whether to only use highly variable genes for PCA
    n_neighbors                    // integer: Number of nearest neighbors to compute
    neighbors_n_pcs                // integer: Number of PCs to use for nearest neighbors
    umap_min_dist                  //   float: Minimum distance between embedded points
    umap_spread                    //   float: Scale of embedded points
    cluster_resolution             //   float: Spot clustering resolution
    rank_genes_method              //  string: Method to use for differential expression testing
    spatial_coord_type             //  string: Type of spatial coordinate system
    spatial_n_neighbors            // integer: Number of spatial neighborhoods
    svg_autocorr_method            // string : Autocorrelation method
    n_top_svgs                     // integer: Number of variable genes to plot
    skip_integration               // boolean: Whether to integrate sdata or not
    integration_method             //  string: Integration method to use
    integration_cluster_resolution //   float: Integration cluster resolution
    multiqc_config                 //    file: /path/to/multiqc/config
    multiqc_logo                   //    file: /path/to/multiqc/logo
    multiqc_methods_description    //    file: /path/to/multiqc/description
    outdir                         //     dir: /path/to/output/directory

    main:
    //
    // WORKFLOW: Run pipeline
    //
    SPATIALVI (
        samplesheet,
        spaceranger_reference,
        spaceranger_probeset,
        hd_bin_size,
        skip_downstream,
        qc_min_counts,
        qc_min_genes,
        qc_min_spots,
        qc_mito_threshold,
        qc_ribo_threshold,
        qc_hb_threshold,
        normalize_target_sum,
        n_highly_variable_genes,
        hvg_flavor,
        n_principal_components,
        pca_use_highly_variable,
        n_neighbors,
        neighbors_n_pcs,
        umap_min_dist,
        umap_spread,
        cluster_resolution,
        rank_genes_method,
        spatial_coord_type,
        spatial_n_neighbors,
        svg_autocorr_method,
        n_top_svgs,
        skip_integration,
        integration_method,
        integration_cluster_resolution,
        multiqc_config,
        multiqc_logo,
        multiqc_methods_description,
        outdir,
    )

    emit:
    multiqc_report = SPATIALVI.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        args,
        params.outdir,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_SPATIALVI (
        params.input,
        params.spaceranger_reference,
        params.spaceranger_probeset,
        params.hd_bin_size,
        params.skip_downstream,
        params.qc_min_counts,
        params.qc_min_genes,
        params.qc_min_spots,
        params.qc_mito_threshold,
        params.qc_ribo_threshold,
        params.qc_hb_threshold,
        params.normalize_target_sum,
        params.n_highly_variable_genes,
        params.hvg_flavor,
        params.n_principal_components,
        params.pca_use_highly_variable,
        params.n_neighbors,
        params.neighbors_n_pcs,
        params.umap_min_dist,
        params.umap_spread,
        params.cluster_resolution,
        params.rank_genes_method,
        params.spatial_coord_type,
        params.spatial_n_neighbors,
        params.svg_autocorr_method,
        params.n_top_svgs,
        params.skip_integration,
        params.integration_method,
        params.integration_cluster_resolution,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir,
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_SPATIALVI.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
