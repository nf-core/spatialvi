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

include { SPATIALVI               } from './workflows/spatialvi'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_spatialvi_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_spatialvi_pipeline'
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
    samplesheet                    // file   : samplesheet read in from --input
    spaceranger_reference          // dir    : /path/to/reference
    spaceranger_probeset           // file   : /path/to/csv
    qc_min_counts                  // integer: Minimum UMIs per spot
    qc_min_genes                   // integer: Minimum genes per spot
    qc_min_spots                   // integer: Minimum spots per gene
    qc_mito_threshold              // float  : Maximum mito. content per spot
    qc_ribo_threshold              // float  : Minimum ribo. content per spot
    qc_hb_threshold                // float  : Maximum haem. content per spot
    cluster_n_hvgs                 // integer: Number of HVGs to use
    cluster_resolution             // float  : Spot clustering resolution
    svg_autocorr_method            // string : Autocorrelation method
    n_top_svgs                     // integer: Number of variable genes to plot
    merge_sdata                    // boolean: Whether to merge sdata or not
    integrate_sdata                // boolean: Whether to integrate sdata or not
    integration_cluster_resolution // float  : Integration cluster resolution
    integration_n_hvgs             // integer: Number of HVGs to integrate with
    multiqc_config                 // file   : /path/to/multiqc/config
    multiqc_logo                   // file   : /path/to/multiqc/logo
    multiqc_methods_description    // file   : /path/to/multiqc/description
    outdir                         // dir    : /path/to/output/directory

    main:
    //
    // WORKFLOW: Run pipeline
    //
    SPATIALVI (
        samplesheet,
        spaceranger_reference,
        spaceranger_probeset,
        qc_min_counts,
        qc_min_genes,
        qc_min_spots,
        qc_mito_threshold,
        qc_ribo_threshold,
        qc_hb_threshold,
        cluster_n_hvgs,
        cluster_resolution,
        svg_autocorr_method,
        n_top_svgs,
        merge_sdata,
        integrate_sdata,
        integration_cluster_resolution,
        integration_n_hvgs,
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
        params.monochrome_logs,
        args,
        params.outdir,
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_SPATIALVI (
        params.input,
        params.spaceranger_reference,
        params.spaceranger_probeset,
        params.qc_min_counts,
        params.qc_min_genes,
        params.qc_min_spots,
        params.qc_mito_threshold,
        params.qc_ribo_threshold,
        params.qc_hb_threshold,
        params.cluster_n_hvgs,
        params.cluster_resolution,
        params.svg_autocorr_method,
        params.n_top_svgs,
        params.merge_sdata,
        params.integrate_sdata,
        params.integration_cluster_resolution,
        params.integration_n_hvgs,
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
