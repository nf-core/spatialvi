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
    // Reference files
    //
    ch_reference = Channel.empty()
    if (params.spaceranger_reference) {
        ch_reference = Channel
            .fromPath ( params.spaceranger_reference, type: "dir", checkIfExists: true )
    } else {
        address = "https://cf.10xgenomics.com/supp/spatial-exp/refdata-gex-mm10-2020-A.tar.gz"
        ch_reference = DOWNLOAD_REFERENCE ( address ).reference
    }

    //
    // Optional: probe set
    //
    ch_probeset = Channel.empty()
    if (params.spaceranger_probeset) {
        ch_probeset = Channel
            .fromPath ( params.spaceranger_probeset, checkIfExists: true )
    } else {
        ch_probeset = file ( 'EMPTY' )
    }

    //
    // Optional: manual alignment file
    //
    ch_manual_alignment = Channel.empty()
    if (params.spaceranger_manual_alignment) {
        ch_manual_alignment = Channel
            .fromPath ( params.spaceranger_manual_alignment, checkIfExists: true )
    } else {
        ch_manual_alignment = file ( 'EMPTY' )
    }

    //
    // Run SpaceRanger count
    //
    SPACERANGER_COUNT (
        ch_input,
        ch_reference,
        ch_probeset,
        ch_manual_alignment
    )
    ch_versions = ch_versions.mix(SPACERANGER_COUNT.out.versions.first())

    emit:
    sr_dir   = SPACERANGER_COUNT.out.sr_dir // channel: [ dir ]
    sr_out   = SPACERANGER_COUNT.out.sr_out // channel: [ val(meta), positions, tissue_lowres_image, tissue_hires_image, scale_factors, barcodes, features, matrix ]

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
