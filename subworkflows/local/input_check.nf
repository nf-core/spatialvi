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
            // collapse technical replicates
            .map { row -> [row.sample, row]}
            .groupTuple()
            .map { id, sample_info -> create_spaceranger_channels(sample_info) }
    } else {
        st_data = SAMPLESHEET_CHECK.out.csv
            .splitCsv ( header: true, sep: ',' )
            .map      { create_visium_channels(it) }
    }

    emit:
    st_data                                   // channel: [ val(meta), [ st data ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

def get_unique(List<LinkedHashMap> sample_info, String key) {
    def val = null
    for (row in sample_info) {
        if (val != null && val != row[key]) {
             exit 1, "ERROR: Please check input samplesheet -> '${key}' is not consistent for all technical replicates of the same sample. Actual: '${row[key]}'. Expected: '${val}'."
        }
        val = row[key]
    }
    return val
}

// Function to get list of [ meta, [ fastq_dir, tissue_hires_image, slide, area ]
def create_spaceranger_channels(List<LinkedHashMap> sample_info) {
    def meta = [:]
    meta.id = get_unique(sample_info, "sample")
    meta.slide = get_unique(sample_info, "slide")
    meta.area = get_unique(sample_info, "area")
    tissue_hires_image = get_unique(sample_info, "tissue_hires_image")
    manual_alignment = get_unique(sample_info, "manual_alignment")
    slidefile = get_unique(sample_info, "slidefile")
    fastq_files = []
    for (row in sample_info) {
        fastq_files.add(file(row.fastq_1, checkIfExists: true))
        fastq_files.add(file(row.fastq_2, checkIfExists: true))
    }

    tissue_hires_image = file(tissue_hires_image, checkIfExists: true)
    manual_alignment = manual_alignment ? file(manual_alignment, checkIfExists: true) : []
    slidefile = slidefile ? file(slidefile, checkIfExists: true) : []
    return [meta, fastq_files, tissue_hires_image, manual_alignment, slidefile]



    // files_to_check = [
    //     "tissue_hires_image": tissue_hires_image,
    //     "tissue_hires_image"
    // ]
    // def raw_meta = []
    // for (entry in row) {
    //     if (entry.key in files_to_check && !file(entry.value).exists()) {
    //         exit 1, "ERROR: Please check input samplesheet -> ${entry.key} file does not exist!\n${entry.value}"
    //     }
    // }
    // if ( row.manual_alignment.isEmpty() ) {
    //     manual_alignment = []
    // } else {
    //     if (!file(row.manual_alignment).exists()) {
    //         exit 1, "ERROR: Please check input samplesheet -> manual_alignment file does not exist!\n${row.manual_alignment}"
    //     }
    //     manual_alignment = file ( row.manual_alignment )
    // }

    // raw_meta = [
    //     meta,
    //     file(row.fastq_dir),
    //     file(row.tissue_hires_image),
    //     manual_alignment
    // ]
    // return raw_meta
}

// Function to get list of [ meta, [ tissue_positions_list, tissue_hires_image, \
// scale_factors, barcodes, features, matrix ] ]
def create_visium_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id = row.sample

    files_to_check = [
        "tissue_positions_list",
        "tissue_lowres_image",
        "tissue_hires_image",
        "scale_factors",
        "barcodes",
        "features",
        "matrix"
    ]
    def processed_meta = []
    for (entry in row) {
        if (entry.key in files_to_check && !file(entry.value).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> ${entry.key} file does not exist!\n${entry.value}"
        }
    }

    processed_meta = [
        meta,
        file(row.tissue_positions_list),
        file(row.tissue_lowres_image),
        file(row.tissue_hires_image),
        file(row.scale_factors),
        file(row.barcodes),
        file(row.features),
        file(row.matrix)
    ]
    return processed_meta
}
