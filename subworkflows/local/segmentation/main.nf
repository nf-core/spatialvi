//
// Cell segmentation of the tissue image
//

include { IMAGE_DOWNSCALE } from "../../../modules/local/image/downscale"
include { CELLPOSE        } from "../../../modules/nf-core/cellpose"
include { STARDIST        } from "../../../modules/nf-core/stardist"
include { MASK_TO_GEOJSON } from "../../../modules/local/segmentation/mask_to_geojson"

workflow SEGMENTATION {

    take:
    ch_image              // channel: [ meta, image ]
    segmentation_method   //  string: cellpose or stardist
    segmentation_downscale // integer: Downscale factor for the image

    main:

    IMAGE_DOWNSCALE (
        ch_image,
        segmentation_downscale
    )

    ch_mask = channel.empty()
    object_class = 'nucleus'
    if (segmentation_method == 'cellpose') {
        CELLPOSE (
            IMAGE_DOWNSCALE.out.image,
            []
        )
        ch_mask = CELLPOSE.out.mask
        object_class = 'cell'
    } else if (segmentation_method == 'stardist') {
        STARDIST (
            IMAGE_DOWNSCALE.out.image,
            ['2D_versatile_he', []]
        )
        ch_mask = STARDIST.out.mask
    } else {
        error("Unsupported segmentation_method '${segmentation_method}': expected 'cellpose' or 'stardist'")
    }

    MASK_TO_GEOJSON (
        ch_mask,
        segmentation_downscale,
        object_class
    )

    emit:
    mask    = ch_mask                    // channel: [ meta, tif ]
    geojson = MASK_TO_GEOJSON.out.geojson // channel: [ meta, geojson ]
}
