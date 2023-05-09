//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {

    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
    if ( params.run_spaceranger ) {
        st_data = SAMPLESHEET_CHECK.out.csv
            .splitCsv ( header: true, sep: ',' )
            .map      { create_spaceranger_channels(it) }
    } else {
        st_data = SAMPLESHEET_CHECK.out.csv
            .splitCsv ( header: true, sep: ',' )
            .map      { create_visium_channels(it) }
    }

    emit:
    st_data                                   // channel: [ val(meta), [ st data ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_dir, tissue_hires_image, slide, area ]
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
        row.area,
    ]
    return array
}

// Function to get list of [ meta, [ tissue_positions_list, tissue_hires_image, \
// scale_factors, barcodes, features, matrix ] ]
def create_visium_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id = row.sample

    def array = []
    if (!file(row.tissue_positions_list).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> tissue_positions_list file does not exist!\n${row.fastq_1}"
    }
    if (!file(row.tissue_lowres_image).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> tissue_lowres_image file does not exist!\n${row.fastq_1}"
    }
    if (!file(row.tissue_hires_image).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> tissue_hires_image file does not exist!\n${row.fastq_1}"
    }
    if (!file(row.scale_factors).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> scale_factors file does not exist!\n${row.fastq_1}"
    }
    if (!file(row.barcodes).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> barcodes file does not exist!\n${row.fastq_1}"
    }
    if (!file(row.features).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> features file does not exist!\n${row.fastq_1}"
    }
    if (!file(row.matrix).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> matrix file does not exist!\n${row.fastq_1}"
    }
    array = [
        meta,
        file(row.tissue_positions_list),
        file(row.tissue_lowres_image),
        file(row.tissue_hires_image),
        file(row.scale_factors),
        file(row.barcodes),
        file(row.features),
        file(row.matrix)
    ]
    return array
}
