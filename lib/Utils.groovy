//
// This file holds several Groovy functions that could be useful for any Nextflow pipeline
//

class Utils {

    public static List DOWNSTREAM_REQUIRED_SPACERANGER_FILES = [
        "raw_feature_bc_matrix.h5",
        "tissue_positions.csv",
        "scalefactors_json.json",
        "tissue_hires_image.png",
        "tissue_lowres_image.png"
    ]

}
