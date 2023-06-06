//
// Raw data processing with Space Ranger
//

include { SPACERANGER_DOWNLOAD_PROBESET  } from '../../modules/local/spaceranger_download_probeset'
include { SPACERANGER_DOWNLOAD_REFERENCE } from '../../modules/local/spaceranger_download_reference'
include { SPACERANGER_COUNT              } from '../../modules/local/spaceranger_count'

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
        ch_reference = file(params.spaceranger_reference, type: "dir", checkIfExists: true)
    } else {
        address = "https://cf.10xgenomics.com/supp/spatial-exp/refdata-gex-GRCh38-2020-A.tar.gz"
        ch_reference = SPACERANGER_DOWNLOAD_REFERENCE ( address ).reference.collect()
    }

    //
    // Optional: probe set
    //
    ch_probeset = Channel.empty()
    if (params.spaceranger_probeset) {
        ch_probeset = file( params.spaceranger_probeset, checkIfExists: true )
    } else {
        ch_probeset = file ( 'EMPTY_PROBESET' )
    }

    //
    // Optional: manual alignment file
    //
    ch_manual_alignment = Channel.empty()
    if (params.spaceranger_manual_alignment) {
        ch_manual_alignment = Channel
            .fromPath ( params.spaceranger_manual_alignment, checkIfExists: true )
    } else {
        ch_manual_alignment = file ( 'EMPTY_ALIGNMENT' )
    }

    //
    // Run Space Ranger count
    //
    SPACERANGER_COUNT (
        ch_st_data,
        ch_reference,
        ch_probeset,
        ch_manual_alignment
    )
    ch_versions = ch_versions.mix(SPACERANGER_COUNT.out.versions.first())

    emit:
    sr_dir   = SPACERANGER_COUNT.out.sr_dir // channel: [ dir ]
    sr_out   = SPACERANGER_COUNT.out.sr_out // channel: [ val(meta), positions, tissue_lowres_image, tissue_hires_image, scale_factors, barcodes, features, matrix ]

    versions = ch_versions                  // channel: [ versions.yml ]
}
