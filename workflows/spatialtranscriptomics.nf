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

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
================================================================================
*/

//
// MODULE: Loaded from modules/local/
//
include { READ_ST_AND_SC_DATA } from '../modules/local/read_st_and_sc_data'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK        } from '../subworkflows/local/input_check'
include { PREPROCESS_SC_DATA } from '../subworkflows/local/preprocess_sc_data'
include { PREPROCESS_ST_DATA } from '../subworkflows/local/preprocess_st_data'
include { SPACERANGER        } from '../subworkflows/local/spaceranger'
include { ST_POSTPROCESSING  } from '../subworkflows/local/st_postprocessing'

/*
================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

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
    // MODULE: Read ST and SC data and save as `anndata`
    //
    READ_ST_AND_SC_DATA (
        ch_st_data
    )
    ch_versions = ch_versions.mix(READ_ST_AND_SC_DATA.out.versions)

    // TODO: Add file manifest or other non-hard-coded path
    //
    // Channel for mitochondrial data
    //
    ch_mito_data = Channel
        .fromPath("ftp://ftp.broadinstitute.org/distribution/metabolic/papers/Pagliarini/MitoCarta2.0/Human.MitoCarta2.0.txt")

    //
    // SUBWORKFLOW: Pre-processing of ST  data
    //
    PREPROCESS_ST_DATA (
        READ_ST_AND_SC_DATA.out.st_counts,
        READ_ST_AND_SC_DATA.out.st_raw,
        ch_mito_data
    )
    ch_versions = ch_versions.mix(PREPROCESS_ST_DATA.out.versions)

    //
    // SUBWORKFLOW (optional): Pre-processing of SC data
    //
    if ( params.single_cell ) {
        PREPROCESS_SC_DATA (
            READ_ST_AND_SC_DATA.out.sc_counts,
            READ_ST_AND_SC_DATA.out.sc_raw,
            ch_mito_data
        )
        ch_versions = ch_versions.mix(PREPROCESS_SC_DATA.out.versions)
    }

    //
    // SUBWORKFLOW: Post-processing and reporting
    //
    ST_POSTPROCESSING (
        PREPROCESS_ST_DATA.out.st_data_norm
    )
    ch_versions = ch_versions.mix(ST_POSTPROCESSING.out.versions)
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
}

/*
================================================================================
    THE END
================================================================================
*/
