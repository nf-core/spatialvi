//
// Raw data processing with Space Ranger
//

include { UNTAR as SPACERANGER_UNTAR_REFERENCE } from "../../../modules/nf-core/untar"
include { SPACERANGER_COUNT                    } from '../../../modules/nf-core/spaceranger/count'

workflow SPACERANGER {

    take:
    ch_data               // channel  : [ val(meta), [ raw st data ] ]
    spaceranger_reference // directory: /path/to/reference
    spaceranger_probeset  // file     : /path/to/csv

    main:

    ch_versions = channel.empty()

    //
    // Reference files
    //
    ch_reference = channel.empty()
    if (spaceranger_reference ==~ /.*\.tar\.gz$/) {
        ref_file = file(spaceranger_reference)
        SPACERANGER_UNTAR_REFERENCE ([
            [id: "reference"],
            ref_file
        ])
        ch_reference = SPACERANGER_UNTAR_REFERENCE.out.untar.map({_meta, ref -> ref})
        ch_versions = ch_versions.mix(SPACERANGER_UNTAR_REFERENCE.out.versions)
    } else {
        ch_reference = file ( spaceranger_reference, type: "dir", checkIfExists: true )
    }

    //
    // Optional: probe set
    //
    ch_probeset = channel.empty()
    if (spaceranger_probeset) {
        ch_probeset = file ( spaceranger_probeset, checkIfExists: true )
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

    emit:
    sr_dir   = SPACERANGER_COUNT.out.outs // channel: [ meta, dir ]
    versions = ch_versions                // channel: [ versions.yml ]
}
