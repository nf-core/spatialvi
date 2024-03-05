//
// Check input samplesheet and get read channels
//

include { UNTAR as UNTAR_SPACERANGER_INPUT } from "../../modules/nf-core/untar"
include { UNTAR as UNTAR_DOWNSTREAM_INPUT  } from "../../modules/nf-core/untar"

workflow INPUT_CHECK {

    take:
    samplesheet // file: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()

    ch_st = Channel.fromPath(samplesheet)
        .splitCsv ( header: true, sep: ',')
        .branch   {
            spaceranger: !it.containsKey("spaceranger_dir")
            downstream: it.containsKey("spaceranger_dir")
        }

    // Space Ranger analysis: --------------------------------------------------

    // Split channel into tarballed and directory inputs
    ch_spaceranger = ch_st.spaceranger
        .map { it -> [it, it.fastq_dir]}
        .branch {
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
    ch_spaceranger_input = ch_spaceranger_combined.map { create_channel_spaceranger(it) }

    // Downstream analysis: ----------------------------------------------------

    // Split channel into tarballed and directory inputs
    ch_downstream = ch_st.downstream
        .map    { create_channel_downstream_tar(it) }
        .branch {
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
    ch_downstream_input = ch_downstream_combined.map { create_channel_downstream(it) }

    emit:
    ch_spaceranger_input   // channel: [ val(meta), [ st data ] ]
    ch_downstream_input    // channel: [ val(meta), [ st data ] ]
    versions = ch_versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ spaceranger_dir ]]
def create_channel_downstream_tar(LinkedHashMap meta) {
    meta['id'] = meta.remove('sample')
    spaceranger_dir = meta.remove('spaceranger_dir')
    return [meta, spaceranger_dir]
}


// Function to get list of [ meta, [ raw_feature_bc_matrix, tissue_positions,
//                                   scalefactors, hires_image, lowres_image ]]
def create_channel_downstream(LinkedHashMap meta) {
    meta["id"] = meta.remove("sample")
    spaceranger_dir = file("${meta.remove('spaceranger_dir')}/**")
    DOWNSTREAM_REQUIRED_SPACERANGER_FILES = [
        "raw_feature_bc_matrix.h5",
        "tissue_positions.csv",
        "scalefactors_json.json",
        "tissue_hires_image.png",
        "tissue_lowres_image.png"
    ]
    for (f in DOWNSTREAM_REQUIRED_SPACERANGER_FILES) {
        if(!spaceranger_dir*.name.contains(f)) {
            error "The specified spaceranger output directory doesn't contain the required file `${f}` for sample `${meta.id}`"
        }
    }
    return [meta, spaceranger_dir]
}

// Function to get list of [ meta, [ fastq_dir, tissue_hires_image, slide, area ]]
def create_channel_spaceranger(LinkedHashMap meta) {
    meta["id"] = meta.remove("sample")

    // Convert a path in `meta` to a file object and return it. If `key` is not contained in `meta`
    // return an empty list which is recognized as 'no file' by nextflow.
    def get_file_from_meta = {key ->
        v = meta.remove(key);
        return v ? file(v) : []
    }

    fastq_dir = meta.remove("fastq_dir")
    fastq_files = file("${fastq_dir}/${meta['id']}*.fastq.gz")
    manual_alignment = get_file_from_meta("manual_alignment")
    slidefile = get_file_from_meta("slidefile")
    image = get_file_from_meta("image")
    cytaimage = get_file_from_meta("cytaimage")
    colorizedimage = get_file_from_meta("colorizedimage")
    darkimage = get_file_from_meta("darkimage")

    if(!fastq_files.size()) {
        error "No `fastq_dir` specified or no samples found in folder."
    }

    check_optional_files = ["manual_alignment", "slidefile", "image", "cytaimage", "colorizedimage", "darkimage"]
    for(k in check_optional_files) {
        if(this.binding[k] && !this.binding[k].exists()) {
            error "File for `${k}` is specified, but does not exist: ${this.binding[k]}."
        }
    }
    if(!(image || cytaimage || colorizedimage || darkimage)) {
        error "Need to specify at least one of 'image', 'cytaimage', 'colorizedimage', or 'darkimage' in the samplesheet"
    }

    return [meta, fastq_files, image, cytaimage, darkimage, colorizedimage, manual_alignment, slidefile]
}

