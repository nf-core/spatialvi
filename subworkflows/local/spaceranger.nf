//
// Raw data processing with SpaceRanger
//

include { DOWNLOAD_PROBESET  } from '../../modules/local/download_probeset'
include { DOWNLOAD_REFERENCE } from '../../modules/local/download_reference'
include { SPACERANGER_COUNT  } from '../../modules/local/spaceranger_count'

workflow SPACERANGER {

    take:
    samplesheet // file: path/to/samplesheet.csv

    main:

    ch_versions = Channel.empty()

    //
    // Read input samplesheet
    //
    ch_input = Channel
        .fromPath ( samplesheet )
        .splitCsv ( header: true, sep: ',' )
        .map      { create_spaceranger_channels(it) }

    //
    // Use a user-provided reference or a default value
    //
    ch_reference = Channel.empty()
    if (params.spaceranger_reference) {
        ch_reference = Channel
            .fromPath ( params.spaceranger_reference, type: "dir", checkIfExists: true )
    } else {
        address = "https://cf.10xgenomics.com/supp/spatial-exp/refdata-gex-mm10-2020-A.tar.gz"
        ch_reference = DOWNLOAD_REFERENCE ( address ).out.reference
    }

    //
    // Use a user-provided probe set or a default value
    //
    ch_probeset = Channel.empty()
    if (params.spaceranger_probeset) {
        ch_probeset = Channel
            .fromPath ( params.spaceranger_probeset, checkIfExists: true )
    } else {
        address = "https://cf.10xgenomics.com/supp/spatial-exp/probeset/Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv"
        ch_probeset = DOWNLOAD_PROBESET ( address ).out.probeset
    }

    //
    // Run SpaceRanger count
    //
    SPACERANGER_COUNT (
        ch_input,
        ch_reference,
        ch_probeset
    )
    ch_versions = ch_versions.mix(SPACERANGER_COUNT.out.versions.first())

    emit:
    sr_dir   = SPACERANGER_COUNT.out.sr_dir // channel: [ dir ]

    versions = ch_versions                  // channel: [ versions.yml ]
}

def create_spaceranger_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id = row.sample

    def array = []
    if (!file(row.fastq_dir).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> fastq_dir directory does not exist!\n${row.fastq_1}"
    }
    if (!file(row.tissue_hires_image).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> tissue_hires_image file does not exist!\n${row.fastq_1}"
    }
    array = [
        meta,
        file(row.fastq_dir),
        file(row.tissue_hires_image),
        row.slide,
        row.area
    ]
    return array
}
