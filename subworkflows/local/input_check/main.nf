//
// Check input samplesheet and get read channels
//

include { UNTAR as UNTAR_SPACERANGER_INPUT } from "../../../modules/nf-core/untar"
include { UNTAR as UNTAR_DOWNSTREAM_INPUT  } from "../../../modules/nf-core/untar"

workflow INPUT_CHECK {

    take:
    samplesheet         // file:    samplesheet read in from --input
    hd_bin_size         // integer: Bin size for Visium HD
    segmentation_method //  string: Segmentation method, if any

    main:

    ch_st = channel.fromPath(samplesheet)
        .splitCsv ( header: true, sep: ',')
        .map      { it -> check_segmentation(it, segmentation_method) }
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

    // Combine extracted and directory inputs into one channel
    ch_spaceranger_combined = UNTAR_SPACERANGER_INPUT.out.untar
        .mix ( ch_spaceranger.dir.map { meta, dir -> [meta, file(dir)] } )
    // Create final meta map and check input existance
    ch_spaceranger_input = ch_spaceranger_combined
        .map { meta, dir -> create_channel_spaceranger(meta, dir, segmentation_method) }

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

    // Combine extracted and directory inputs into one channel
    ch_downstream_combined = UNTAR_DOWNSTREAM_INPUT.out.untar
        .mix ( ch_downstream.dir )
        .map { meta, dir -> [meta, dir] }

    // Create final meta map and check input file existence
    ch_downstream_input = ch_downstream_combined
        .map { it -> check_downstream_dir(it, hd_bin_size, segmentation_method) }

    // Tissue images from both input modes
    ch_image = ch_spaceranger_input
        .map { it -> [it[0], it[2]] }
        .mix ( ch_downstream_input.map { meta, _dir -> [meta, get_file_from_meta(meta, "image")] } )
        .filter { _meta, image -> image as boolean }

    // Segmentations supplied in the samplesheet, in full-resolution pixels
    ch_segmentation_file = ch_spaceranger_input
        .map { it -> it[0] }
        .mix ( ch_downstream_input.map { meta, _dir -> meta } )
        .map { meta -> [meta, get_file_from_meta(meta, "segmentation_file")] }
        .filter { _meta, segmentation -> segmentation as boolean }

    emit:
    ch_spaceranger_input   // channel: [ val(meta), [ st data ] ]
    ch_downstream_input    // channel: [ val(meta), [ st data ] ]
    ch_image               // channel: [ val(meta), image ]
    ch_segmentation_file   // channel: [ val(meta), segmentation ]
}

// Function: `segmentation_file` and `--segmentation_method` are mutually exclusive
def check_segmentation(meta, segmentation_method) {
    if (segmentation_method && meta.get('segmentation_file')) {
        error "Sample '${meta.sample}' has a `segmentation_file`, which cannot be combined with `--segmentation_method`."
    }
    return meta
}

// Function: normalize meta only (no filesystem checks)
def create_channel_downstream(meta) {
    meta['id'] = meta.get('sample')
    def spaceranger_dir = meta.get('spaceranger_dir')
    return [meta, spaceranger_dir]
}

// Function: validate that required files exist in the dir
def check_downstream_dir(input, hd_bin_size, segmentation_method) {
    def (meta, spaceranger_dir) = input

    // Check for existence of optional files
    def image = get_file_from_meta(meta, "image")
    def segmentation_file = get_file_from_meta(meta, "segmentation_file")
    ['image': image, 'segmentation_file': segmentation_file].each { k, f ->
        if (f && !f.exists()) {
            error "File for `${k}` is specified, but does not exist: ${f}."
        }
    }

    // Non-HD SpaceRanger output required files
    def classic_required_files = [
        "filtered_feature_bc_matrix.h5",
        "tissue_positions.csv",
        "scalefactors_json.json",
        "tissue_hires_image.png",
        "tissue_lowres_image.png"
    ]
    def dir_file_objs = file("${spaceranger_dir}/**")
    def classic_files_present = classic_required_files
        .every { f -> dir_file_objs*.name.contains(f) }

    // Visium HD required files (for specified bin size)
    def hd_bin_str = "square_${String.format('%03d', hd_bin_size)}um"
    def hd_required_files = [
        "feature_slice.h5",
        "binned_outputs/${hd_bin_str}/filtered_feature_bc_matrix.h5",
        "binned_outputs/${hd_bin_str}/spatial/scalefactors_json.json",
        "binned_outputs/${hd_bin_str}/spatial/tissue_hires_image.png",
        "binned_outputs/${hd_bin_str}/spatial/tissue_lowres_image.png",
        "binned_outputs/${hd_bin_str}/spatial/tissue_positions.parquet"
    ]
    def hd_files_present = hd_required_files.every { f -> file("${spaceranger_dir}/${f}").exists() }

    if (!(classic_files_present || hd_files_present)) {
        error "The specified spaceranger output directory for sample '${meta.id}' does not contain all required files for either classic Visium: ${classic_required_files.join(', ')} or Visium HD bin size ${hd_bin_size}: ${hd_required_files.join(', ')}."
    }

    // Route flags: bin2cell groups the bins of processed data into cells
    meta.hd          = hd_files_present
    meta.cells       = (segmentation_file || segmentation_method) as boolean
    meta.bin_to_cell = (segmentation_file || segmentation_method in ['cellpose', 'stardist']) as boolean

    // Cell segmentation is Visium HD only
    if (meta.cells && !hd_files_present) {
        error "Sample '${meta.id}': cell segmentation was requested but '${spaceranger_dir}' is classic Visium output. Segmentation is only supported for Visium HD."
    }
    check_segmentation_image(meta, image, segmentation_method)
    if (meta.bin_to_cell) {
        def bins_2um = "binned_outputs/square_002um/filtered_feature_bc_matrix.h5"
        if (!file("${spaceranger_dir}/${bins_2um}").exists()) {
            error "Sample '${meta.id}': grouping bins into cells needs the 2 µm bins, but '${spaceranger_dir}/${bins_2um}' is missing."
        }
    }
    if (segmentation_method == 'spaceranger' && !file("${spaceranger_dir}/segmented_outputs").isDirectory()) {
        error "Sample '${meta.id}': --segmentation_method spaceranger needs '${spaceranger_dir}/segmented_outputs', which is missing. Space Ranger 4.0 or later produces it for Visium HD."
    }

    return [meta, spaceranger_dir]
}

// Function: cellpose and stardist segment the samplesheet `image`
def check_segmentation_image(meta, image, segmentation_method) {
    if (segmentation_method in ['cellpose', 'stardist'] && !image) {
        error "Sample '${meta.id}': `--segmentation_method ${segmentation_method}` needs an `image` in the samplesheet."
    }
}

// Function to convert a path in `meta` to a file object and return it. If key
// `k` is not contained in `meta` return an empty list which is recognized as
// 'no file' by Nextflow.
def get_file_from_meta(meta, k) {
    def v = meta.get(k)
    return v ? file(v) : []
}

// Function to get [ meta, [fastq_dir, tissue_hires_image, slide, area] ] list
def create_channel_spaceranger(meta, fastq_dir, segmentation_method) {
    meta["id"] = meta.get("sample")
    def slide = meta.get("slide")
    def area = meta.get("area")

    // Resolve symlinks for local filesystem paths only
    def scheme = fastq_dir.toUri().getScheme()
    if (scheme == null || scheme == 'file') {
        fastq_dir = fastq_dir.toRealPath()
    }

    def fastq_files = fastq_dir.listFiles().findAll { file ->
        file.name.endsWith('.fastq.gz')
    }

    def manual_alignment = get_file_from_meta(meta, "manual_alignment")
    def slidefile = get_file_from_meta(meta, "slidefile")
    def image = get_file_from_meta(meta, "image")
    def cytaimage = get_file_from_meta(meta, "cytaimage")
    def colorizedimage = get_file_from_meta(meta, "colorizedimage")
    def darkimage = get_file_from_meta(meta, "darkimage")
    def segmentation_file = get_file_from_meta(meta, "segmentation_file")

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
        'darkimage': darkimage,
        'segmentation_file': segmentation_file
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

    // Route flags: Space Ranger segments raw data itself, so never bin2cell
    meta.hd          = (slide?.startsWith('H1-') || (slidefile && slidefile.name.startsWith('H1-'))) as boolean
    meta.cells       = (segmentation_file || segmentation_method) as boolean
    meta.bin_to_cell = false

    // Cell segmentation is Visium HD only; allowed when HD is undetectable
    if (meta.cells && slide && !slide.startsWith('H1-')) {
        error "Sample '${meta.id}': cell segmentation was requested but slide '${slide}' is classic Visium. Segmentation is only supported for Visium HD."
    }
    check_segmentation_image(meta, image, segmentation_method)

    // `segmentation_file` is emitted on its own channel, not in this tuple
    return [meta, fastq_files, image, slide, area, cytaimage, darkimage, colorizedimage, manual_alignment, slidefile]
}
