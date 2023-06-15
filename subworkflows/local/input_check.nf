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
def create_spaceranger_channels(LinkedHashMap meta) {
    meta["id"] = meta["sample"]
    meta.remove("sample")

    fastq_dir = meta.remove("fastq_dir")

    files_to_check = [
        "fastq_dir",
        "image",
        "cytaimage",
        "colorizedimage",
        "darkimage"
    ]
    def raw_meta = []
    for (entry in row) {
        if (entry.key in files_to_check && !file(entry.value).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> ${entry.key} file does not exist!\n${entry.value}"
        }
    }
    if ( row.manual_alignment.isEmpty() ) {
        manual_alignment = []
    } else {
        if (!file(row.manual_alignment).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> manual_alignment file does not exist!\n${row.manual_alignment}"
        }
        manual_alignment = file ( row.manual_alignment )
    }

    raw_meta = [
        meta,
        file(row.fastq_dir),
        file(row.tissue_hires_image),
        row.slide,
        row.area,
        manual_alignment
    ]
    return raw_meta
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
