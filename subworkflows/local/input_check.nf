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

    def get_file_from_meta = {key ->
         v = meta.remove(key);
         return v ? file(v) : []
    }

    fastq_dir = meta.remove("fastq_dir")
    manual_alignment = get_file_from_meta("manual_alignment")
    println("${fastq_dir}/${meta['id']}*.fastq.gz")
    fastq_files = Channel.fromPath("${fastq_dir}/${meta['id']}*.fastq.gz").view().to_list().view()
    slidefile = get_file_from_meta("slidefile")
    image = get_file_from_meta("image")
    cytaimage = get_file_from_meta("cytaimage")
    colorizedimage = get_file_from_meta("colorizedimage")
    darkimage = get_file_from_meta("darkimage")

    // if(!fastq_files.length) {
    //     error "No `fastq_dir` specified or no samples found in folder."
    // } else {
    //     log.info "${fastq_files.length} FASTQ files found for sample ${meta['id']}."
    // }

    // if(manual_alignment && !manual_alignment.exist()) {
    //     error "Manual alignment file does not exist: ${manual_alignment}"
    // }
    // if(slidefile && !slidefile.exist()) {
    //     error "Slidefile does not exist: ${manual_alignment}"
    // }
    // if(!(image || cytaimage || colorizedimage || darkimage)) {
    //     error "Need to specify at least one of 'image', 'cytaimage', 'colorizedimage', or 'darkimage' in the samplesheet"
    // }
    // println(this.binding)

    return [meta, fastq_files, image, cytaimage, darkimage, colorizedimage, manual_alignment, slidefile]
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
