/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
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
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

include { ST_PREPARE_DATA          } from '../subworkflows/local/stPrepareData'
include { ST_LOAD_PREPROCESS_DATA  } from '../subworkflows/local/stLoadPreprocessData'
include { ST_MISCELLANEOUS_TOOLS   } from '../subworkflows/local/stMiscellaneousTools'
include { ST_POSTPROCESSING        } from '../subworkflows/local/stPostprocessing'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
//include { FASTQC                      } from '../modules/nf-core/modules/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

outdir = "${launchDir}/${params.outdir}"

Channel
.from(file(params.input))
.splitCsv(header:true, sep:',')
.map{ prep_input_csv_files(it) }
.set{sample_ids}

import org.apache.commons.csv.CSVPrinter
import org.apache.commons.csv.CSVFormat
import groovy.json.JsonOutput
import groovy.json.JsonSlurper

def loadFromURL(String sample_id, String data_dir, String pref_dir, List files) {
    for(cfile: files) {    
        def lfile = cfile.split('/')
        if (lfile.size()>1) {
            subdir = String.join('/', lfile[0..lfile.size()-2]) + '/'
        } else {
            subdir = ''
        }
        
        filename = lfile[lfile.size()-1]        
        url = data_dir + subdir + filename
        savedir = outdir + '/' + 'dataset-' + sample_id + '/' + pref_dir + '/' + subdir
    
        File filedir = new File(savedir)
        filedir.mkdirs()
     
        URL urlobj = new URL(url)
        File fileload = new File(savedir + filename)
        if (!fileload.exists()) {
            println 'Loading: ' + url
            fileload.bytes = urlobj.bytes
        }
    }
    
    return outdir + '/' + 'dataset-' + sample_id + '/' + pref_dir + '/'
}

def prep_input_csv_files(LinkedHashMap row) {
    
    files_st = ['spatial/detected_tissue_image.jpg',
                'spatial/scalefactors_json.json',
                'spatial/tissue_hires_image.png',
                'spatial/tissue_positions_list.csv',
                'spatial/aligned_fiducials.jpg',
                'spatial/tissue_lowres_image.png',
                'raw_feature_bc_matrix/features.tsv.gz',
                'raw_feature_bc_matrix/barcodes.tsv.gz',
                'raw_feature_bc_matrix/matrix.mtx.gz']
             
    files_sc = ['features.tsv.gz',
                'barcodes.tsv.gz',
                'matrix.mtx.gz']
                    
    if (row.st_data_dir[0..3]=='http') {
        row.st_data_dir = loadFromURL(row.sample_id, row.st_data_dir, 'ST', files_st)
    }
    
    if (row.sc_data_dir[0..3]=='http') {
        row.sc_data_dir = loadFromURL(row.sample_id, row.sc_data_dir, 'SC', files_sc)
    }
    
    def fileName = String.format("%s/sample_%s", outdir, row.sample_id)    
    def FILE_HEADER = row.keySet() as String[];
    
    new File(fileName + ".csv").withWriter { fileWriter ->
        def csvFilePrinter = new CSVPrinter(fileWriter, CSVFormat.DEFAULT)
        csvFilePrinter.printRecord(FILE_HEADER)
        csvFilePrinter.printRecord(row.values())
    }
       
    File file = new File(fileName + ".json")
    file.write(JsonOutput.toJson(row))
    
    return row.sample_id
}



workflow ST {

    ST_PREPARE_DATA(             sample_ids,                   outdir )
    
    ST_LOAD_PREPROCESS_DATA(     ST_PREPARE_DATA.out,          outdir )
                          
    ST_MISCELLANEOUS_TOOLS(      ST_LOAD_PREPROCESS_DATA.out,  outdir )
    
    ST_POSTPROCESSING(           ST_MISCELLANEOUS_TOOLS.out,   outdir )
    
    //ST_POSTPROCESSING.out.view()
    
}


workflow ST_PROPER {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowSpatialtranscriptomics.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
