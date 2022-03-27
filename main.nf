#!/usr/bin/env nextflow
/*
================================================================================
    nf-core/spatialtranscriptomics
================================================================================
    Github : https://github.com/nf-core/spatialtranscriptomics
    Website: https://nf-co.re/st
    Slack  : https://nfcore.slack.com/channels/st
--------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
================================================================================
    GENOME PARAMETER VALUES
================================================================================
*/

params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

/*
================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
================================================================================
    NAMED WORKFLOW FOR PIPELINE
================================================================================
*/

include { ST } from './workflows/spatialtranscriptomics'

//
// WORKFLOW: Run main nf-core/spatialtranscriptomics analysis pipeline
//
workflow NFCORE_ST {
    ST ()
}

/*
================================================================================
    RUN ALL WORKFLOWS
================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_ST ()
}

/*
================================================================================
    THE END
================================================================================
*/
