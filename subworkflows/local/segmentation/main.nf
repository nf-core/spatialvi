//
// Cell segmentation of the tissue image
//

include { IMAGE_DOWNSCALE } from "../../../modules/local/image/downscale"
include { CELLPOSE        } from "../../../modules/nf-core/cellpose"
include { MASK_TO_GEOJSON } from "../../../modules/local/segmentation/mask_to_geojson"

workflow SEGMENTATION {

    take:
    ch_image              // channel: [ meta, image ]
    segmentation_method   //  string: cellpose
    segmentation_downscale // integer: Downscale factor for the image

    main:

    IMAGE_DOWNSCALE (
        ch_image,
        segmentation_downscale
    )

    if (segmentation_method != 'cellpose') {
        error("Unsupported segmentation_method '${segmentation_method}': expected 'cellpose'")
    }

    CELLPOSE (
        IMAGE_DOWNSCALE.out.image,
        []
    )

    MASK_TO_GEOJSON (
        CELLPOSE.out.mask,
        segmentation_downscale,
        'cell'
    )

    emit:
    mask    = CELLPOSE.out.mask          // channel: [ meta, tif ]
    geojson = MASK_TO_GEOJSON.out.geojson // channel: [ meta, geojson ]
}
