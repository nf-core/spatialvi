/*
================================================================================
    VALIDATE INPUTS
================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowSpatialtranscriptomics.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
log.info """\
         Project directory:  ${projectDir}
         """
         .stripIndent()


def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
================================================================================
    CONFIG FILES
================================================================================
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
================================================================================
*/

//
// MODULE: Loaded from modules/local/
//
include { ST_READ_DATA } from '../modules/local/st_read_data'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK    } from '../subworkflows/local/input_check'
include { SPACERANGER    } from '../subworkflows/local/spaceranger'
include { ST_PREPROCESS  } from '../subworkflows/local/st_preprocess'
include { ST_POSTPROCESS } from '../subworkflows/local/st_postprocess'

/*
================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
================================================================================
    RUN MAIN WORKFLOW
================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

//
// Spatial transcriptomics workflow
//
workflow ST {

    // TODO: Collect versions for all modules/subworkflows
    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // SUBWORKFLOW: SpaceRanger raw data processing
    //
    if ( params.run_spaceranger ) {
        SPACERANGER (
            params.spaceranger_input
        )
        ch_st_data = SPACERANGER.out.sr_out
        ch_versions = ch_versions.mix(SPACERANGER.out.versions)
    } else {
        ch_st_data = INPUT_CHECK.out.reads
    }

    //
    // MODULE: Read ST data and save as `anndata`
    //
    ST_READ_DATA (
        ch_st_data
    )
    ch_versions = ch_versions.mix(ST_READ_DATA.out.versions)

    // TODO: Add file manifest or other non-hard-coded path
    //
    // Channel for mitochondrial data
    //
    ch_mito_data = Channel
        .fromPath("ftp://ftp.broadinstitute.org/distribution/metabolic/papers/Pagliarini/MitoCarta2.0/Human.MitoCarta2.0.txt")

    //
    // SUBWORKFLOW: Pre-processing of ST  data
    //
    ST_PREPROCESS (
        ST_READ_DATA.out.st_raw,
        ch_mito_data
    )
    ch_versions = ch_versions.mix(ST_PREPROCESS.out.versions)

    //
    // SUBWORKFLOW (optional): Pre-processing of SC data
    //
    if ( params.single_cell ) {
        PREPROCESS_SC_DATA (
            ST_READ_DATA.out.sc_counts,
            ST_READ_DATA.out.sc_raw,
            ch_mito_data
        )
        ch_versions = ch_versions.mix(PREPROCESS_SC_DATA.out.versions)
    }

    //
    // SUBWORKFLOW: Post-processing and reporting
    //
    ST_POSTPROCESS (
        ST_PREPROCESS.out.st_data_norm
    )
    ch_versions = ch_versions.mix(ST_POSTPROCESS.out.versions)
}

/*
================================================================================
    COMPLETION EMAIL AND SUMMARY
================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
================================================================================
    THE END
================================================================================
*/
