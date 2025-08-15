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
    ch_downstream_input = ch_downstream_combined.map { it -> create_channel_downstream_tar(it) }

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

// Function to get list of [ meta, [ fastq_dir, tissue_hires_image, slide, area ]]
def create_channel_spaceranger(meta) {
    meta["id"] = meta.remove("sample")
    def slide = meta.remove("slide")
    def area = meta.remove("area")
    def fastq_dir = meta.remove("fastq_dir")
    def fastq_files = file("${fastq_dir}/${meta['id']}*.fastq.gz")

    // Convert a path in `meta` to a file object and return it. If key `k` is
    // not contained in `meta` return an empty list which is recognized as 'no
    // file' by Nextflow.
    def get_file_from_meta = { k ->
        def v = meta.remove(k)
        return v ? file(v) : []
    }
    def manual_alignment = get_file_from_meta("manual_alignment")
    def slidefile = get_file_from_meta("slidefile")
    def image = get_file_from_meta("image")
    def cytaimage = get_file_from_meta("cytaimage")
    def colorizedimage = get_file_from_meta("colorizedimage")
    def darkimage = get_file_from_meta("darkimage")

    if(!fastq_files.size()) {
        error "No `fastq_dir` specified or no samples found in folder."
    }

    // Check for existance of optional files
    def optional_files = [
        'manual_alignment': manual_alignment,
        'slidefile': slidefile,
        'image': image,
        'cytaimage': cytaimage,
        'colorizedimage': colorizedimage,
        'darkimage': darkimage
    ]
    optional_files.each { k, f ->
        if(f && !f.exists()) {
            error "File for `${k}` is specified, but does not exist: ${f}."
        }
    }

    // Check that at least one type of image is specified
    if(!(image || cytaimage || colorizedimage || darkimage)) {
        error "Need to specify at least one of 'image', 'cytaimage', 'colorizedimage', or 'darkimage' in the samplesheet"
    }

    return [meta, fastq_files, image, slide, area, cytaimage, darkimage, colorizedimage, manual_alignment, slidefile]
}
