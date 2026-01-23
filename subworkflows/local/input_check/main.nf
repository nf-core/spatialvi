//
// Check input samplesheet and get read channels
//

include { UNTAR as UNTAR_SPACERANGER_INPUT } from "../../../modules/nf-core/untar"
include { UNTAR as UNTAR_DOWNSTREAM_INPUT  } from "../../../modules/nf-core/untar"

workflow INPUT_CHECK {

    take:
    samplesheet // file:    samplesheet read in from --input
    hd_bin_size // integer: Bin size for Visium HD

    main:

    ch_versions = channel.empty()

    ch_st = channel.fromPath(samplesheet)
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
            tar: it[1] ==~ /.*\.tar(\.gz)?$/
            dir: true
        }

    // Extract tarballed inputs
    UNTAR_SPACERANGER_INPUT ( ch_spaceranger.tar )
    ch_versions = ch_versions.mix(UNTAR_SPACERANGER_INPUT.out.versions)

    // Combine extracted and directory inputs into one channel
    ch_spaceranger_combined = UNTAR_SPACERANGER_INPUT.out.untar
        .mix ( ch_spaceranger.dir.map { meta, dir -> [meta, file(dir)] } )
    // Create final meta map and check input existance
    ch_spaceranger_input = ch_spaceranger_combined.map { meta, dir -> create_channel_spaceranger(meta, dir) }

    // Downstream analysis: ----------------------------------------------------

    // Split channel into tarballed and directory inputs
    ch_downstream = ch_st.downstream
        .map    { it -> create_channel_downstream(it) }
        .branch { it ->
            tar: it[1] ==~ /.*\.tar(\.gz)?$/
            dir: true
        }

    // Extract tarballed inputs
    UNTAR_DOWNSTREAM_INPUT ( ch_downstream.tar )
    ch_versions = ch_versions.mix(UNTAR_DOWNSTREAM_INPUT.out.versions)

    // Combine extracted and directory inputs into one channel
    ch_downstream_combined = UNTAR_DOWNSTREAM_INPUT.out.untar
        .mix ( ch_downstream.dir )
        .map { meta, dir -> [meta, dir] }

    // Create final meta map and check input file existence
    ch_downstream_input = ch_downstream_combined.map { it -> check_downstream_dir(it, hd_bin_size) }

    emit:
    ch_spaceranger_input   // channel: [ val(meta), [ st data ] ]
    ch_downstream_input    // channel: [ val(meta), [ st data ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
// Function: normalize meta only (no filesystem checks)
def create_channel_downstream(meta) {
    meta['id'] = meta.get('sample')
    def spaceranger_dir = meta.get('spaceranger_dir')
    return [meta, spaceranger_dir]
}

// Function: validate that required files exist in the dir
def check_downstream_dir(input, hd_bin_size) {
    def (meta, spaceranger_dir) = input

    // Non-HD SpaceRanger output required files
    def classic_required_files = [
        "raw_feature_bc_matrix.h5",
        "tissue_positions.csv",
        "scalefactors_json.json",
        "tissue_hires_image.png",
        "tissue_lowres_image.png"
    ]
    def dir_file_objs = file("${spaceranger_dir}/**")
    def classic_files_present = classic_required_files.every { f -> dir_file_objs*.name.contains(f) }

    // Visium HD binned output required files (for specified bin size)
    def hd_required_files = [
        "raw_feature_bc_matrix.h5",
        "spatial/scalefactors_json.json",
        "spatial/tissue_hires_image.png",
        "spatial/tissue_lowres_image.png",
        "spatial/tissue_positions.parquet"
    ]
    def hd_dir = file("${spaceranger_dir}/binned_outputs/square_${String.format('%03d', hd_bin_size)}um")
    def hd_files_present = hd_required_files.every { f -> file("${hd_dir}/${f}").exists() }

    if (!(classic_files_present || hd_files_present)) {
        error "The specified spaceranger output directory for sample '${meta.id}' does not contain all required files for either classic Visium: ${classic_required_files.join(', ')} or Visium HD bin size ${hd_bin_size}: ${hd_required_files.join(', ')}."
    }

    return [meta, spaceranger_dir]
}

// Function to get list of [ meta, [ fastq_dir, tissue_hires_image, slide, area ]]
def create_channel_spaceranger(meta, fastq_dir) {
    meta["id"] = meta.get("sample")
    def slide = meta.get("slide")
    def area = meta.get("area")

    // Resolve symlinks for local filesystem paths only
    def scheme = fastq_dir.toUri().getScheme()
    if (scheme == null || scheme == 'file') {
        fastq_dir = fastq_dir.toRealPath() // resolve symlink (if applicable)
    }

    def fastq_files = fastq_dir.listFiles().findAll { file ->
        file.name.endsWith('.fastq.gz')
    }

    // Convert a path in `meta` to a file object and return it. If key `k` is
    // not contained in `meta` return an empty list which is recognized as 'no
    // file' by Nextflow.
    def get_file_from_meta = { k ->
        def v = meta.get(k)
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
