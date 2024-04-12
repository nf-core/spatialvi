//
// Raw data processing with Space Ranger
//

include { UNTAR as SPACERANGER_UNTAR_REFERENCE } from "../../modules/nf-core/untar"
include { SPACERANGER_COUNT                    } from '../../modules/nf-core/spaceranger/count'

workflow SPACERANGER {

    take:
    ch_data // channel: [ val(meta), [ raw st data ] ]

    main:

    ch_versions = Channel.empty()

    //
    // Reference files
    //
    ch_reference = Channel.empty()
    if (params.spaceranger_reference ==~ /.*\.tar\.gz$/) {
        ref_file = file(params.spaceranger_reference)
        SPACERANGER_UNTAR_REFERENCE ([
            [id: "reference"],
            ref_file
        ])
        ch_reference = SPACERANGER_UNTAR_REFERENCE.out.untar.map({meta, ref -> ref})
        ch_versions = ch_versions.mix(SPACERANGER_UNTAR_REFERENCE.out.versions)
    } else {
        ch_reference = file ( params.spaceranger_reference, type: "dir", checkIfExists: true )
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
        ch_data,
        ch_reference,
        ch_probeset
    )
    ch_versions = ch_versions.mix(SPACERANGER_COUNT.out.versions.first())

    emit:
    sr_dir   = SPACERANGER_COUNT.out.outs
    versions = ch_versions                  // channel: [ versions.yml ]
}
