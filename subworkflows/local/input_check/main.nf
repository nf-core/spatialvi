//
// Check input samplesheet and get read channels
//

include { UNTAR as UNTAR_SPACERANGER_INPUT } from "../../../modules/nf-core/untar"
include { UNTAR as UNTAR_DOWNSTREAM_INPUT  } from "../../../modules/nf-core/untar"

workflow INPUT_CHECK {

    take:
    samplesheet // file: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()

    ch_st = Channel.fromPath(samplesheet)
        .splitCsv ( header: true, sep: ',')
        .branch   { it ->
            spaceranger: !it.containsKey("spaceranger_dir")
            downstream: it.containsKey("spaceranger_dir")
        }

    // Space Ranger analysis: --------------------------------------------------

    // Split channel into tarballed and directory inputs
    ch_spaceranger = ch_st.spaceranger
        .map { it -> [it, it.fastq_dir]}
        .branch { it ->
            tar: it[1].contains(".tar.gz")
            dir: !it[1].contains(".tar.gz")
        }

    // Extract tarballed inputs
    UNTAR_SPACERANGER_INPUT ( ch_spaceranger.tar )
    ch_versions = ch_versions.mix(UNTAR_SPACERANGER_INPUT.out.versions)

    // Combine extracted and directory inputs into one channel
    ch_spaceranger_combined = UNTAR_SPACERANGER_INPUT.out.untar
        .mix ( ch_spaceranger.dir )
        .map { meta, dir -> meta + [fastq_dir: dir] }

    // Create final meta map and check input existance
    ch_spaceranger_input = ch_spaceranger_combined.map { it -> create_channel_spaceranger(it) }

    // Downstream analysis: ----------------------------------------------------

    // Split channel into tarballed and directory inputs
    ch_downstream = ch_st.downstream
        .map    { it -> create_channel_downstream_tar(it) }
        .branch { it ->
            tar: it[1].contains(".tar.gz")
            dir: !it[1].contains(".tar.gz")
        }

    // Extract tarballed inputs
    UNTAR_DOWNSTREAM_INPUT ( ch_downstream.tar )
    ch_versions = ch_versions.mix(UNTAR_DOWNSTREAM_INPUT.out.versions)

    // Combine extracted and directory inputs into one channel
    ch_downstream_combined = UNTAR_DOWNSTREAM_INPUT.out.untar
        .mix ( ch_downstream.dir )
        .map { meta, dir -> [sample: meta.id, spaceranger_dir: dir] }

    // Create final meta map and check input file existance
    ch_downstream_input = ch_downstream_combined.map { it -> create_channel_downstream(it) }

    emit:
    ch_spaceranger_input   // channel: [ val(meta), [ st data ] ]
    ch_downstream_input    // channel: [ val(meta), [ st data ] ]
    versions = ch_versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ spaceranger_dir ]]
def create_channel_downstream_tar(meta) {
    meta['id'] = meta.remove('sample')
    def spaceranger_dir = meta.remove('spaceranger_dir')
    return [meta, spaceranger_dir]
}


// Function to get list of [ meta, [ raw_feature_bc_matrix, tissue_positions,
//                                   scalefactors, hires_image, lowres_image ]]
def create_channel_downstream(meta) {
    meta["id"] = meta.remove("sample")
    def spaceranger_dir = file("${meta.remove('spaceranger_dir')}/**")
    def DOWNSTREAM_REQUIRED_SPACERANGER_FILES = [
        "raw_feature_bc_matrix.h5",
        "tissue_positions.csv",
        "scalefactors_json.json",
        "tissue_hires_image.png",
        "tissue_lowres_image.png"
    ]
    DOWNSTREAM_REQUIRED_SPACERANGER_FILES.each { f ->
        if(!spaceranger_dir*.name.contains(f)) {
            error "The specified spaceranger output directory doesn't contain the required file `${f}` for sample `${meta.id}`"
        }
    }
    return [meta, spaceranger_dir]
}

// Function to get list of [ meta, [ fastq_dir, tissue_hires_image, slide, area ]]
def create_channel_spaceranger(meta) {
    meta["id"] = meta.remove("sample")
    def slide = meta.remove("slide")
    def area = meta.remove("area")

    // Convert a path in `meta` to a file object and return it. If `key` is not contained in `meta`
    // return an empty list which is recognized as 'no file' by nextflow.
    def get_file_from_meta = {key ->
        def v = meta.remove(key);
        return v ? file(v) : []
    }

    def fastq_dir = meta.remove("fastq_dir")
    def fastq_files = file("${fastq_dir}/${meta['id']}*.fastq.gz")
    def manual_alignment = get_file_from_meta("manual_alignment")
    def slidefile = get_file_from_meta("slidefile")
    def image = get_file_from_meta("image")
    def cytaimage = get_file_from_meta("cytaimage")
    def colorizedimage = get_file_from_meta("colorizedimage")
    def darkimage = get_file_from_meta("darkimage")

    if(!fastq_files.size()) {
        error "No `fastq_dir` specified or no samples found in folder."
    }

    def check_optional_files = ["manual_alignment", "slidefile", "image", "cytaimage", "colorizedimage", "darkimage"]
    check_optional_files.each { k ->
        if(this.binding[k] && !this.binding[k].exists()) {
            error "File for `${k}` is specified, but does not exist: ${this.binding[k]}."
        }
    }
    if(!(image || cytaimage || colorizedimage || darkimage)) {
        error "Need to specify at least one of 'image', 'cytaimage', 'colorizedimage', or 'darkimage' in the samplesheet"
    }

    return [meta, fastq_files, image, slide, area, cytaimage, darkimage, colorizedimage, manual_alignment, slidefile]
}
