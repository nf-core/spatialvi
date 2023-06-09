//
// Raw data processing with Space Ranger
//

include { UNTAR as SPACERANGER_DOWNLOAD_REFERENCE } from "../../modules/nf-core/untar"
include { SPACERANGER_COUNT              } from '../../modules/nf-core/spaceranger/count'

workflow SPACERANGER {

    take:
    ch_st_data // channel: [ val(meta), [ raw st data ] ]

    main:

    ch_versions = Channel.empty()

    //
    // Reference files
    //
    ch_reference = Channel.empty()
    if (params.spaceranger_reference) {
        ch_reference = file ( params.spaceranger_reference, type: "dir", checkIfExists: true )
    } else {
        SPACERANGER_DOWNLOAD_REFERENCE ([
            [id: "refdata-gex-GRCh38-2020-A"],
            file("https://cf.10xgenomics.com/supp/spatial-exp/refdata-gex-GRCh38-2020-A.tar.gz")
        ])
        ch_reference = SPACERANGER_DOWNLOAD_REFERENCE.out.untar.map({meta, ref -> ref})
        ch_versions = ch_versions.mix(SPACERANGER_DOWNLOAD_REFERENCE.out.versions)
    }

    //
    // Optional: probe set
    //
    ch_probeset = Channel.empty()
    if (params.spaceranger_probeset) {
        ch_probeset = file ( params.spaceranger_probeset, checkIfExists: true )
    } else {
        ch_probeset = []
    }

    //
    // Run Space Ranger count
    //
    SPACERANGER_COUNT (
        ch_st_data,
        ch_reference,
        ch_probeset
    )
    ch_versions = ch_versions.mix(SPACERANGER_COUNT.out.versions.first())

    emit:
    sr_dir   = SPACERANGER_COUNT.out.sr_dir // channel: [ dir ]
    sr_out   = SPACERANGER_COUNT.out.sr_out // channel: [ val(meta), positions, tissue_lowres_image, tissue_hires_image, scale_factors, barcodes, features, matrix ]

    versions = ch_versions                  // channel: [ versions.yml ]
}
