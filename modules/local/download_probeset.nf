//
// Download SpaceRanger probeset
//

process DOWNLOAD_PROBESET {

    tag "${name}"
    label "process_low"

    conda (params.enable_conda ? "conda-forge::gnu-wget=1.18" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--hed695b0_4' :
        'quay.io/biocontainers/gnu-wget:1.18--hed695b0_4' }"

    input:
    val(address)

    output:
    path("*.csv"), emit: probeset

    script:
    probeset = address.tokenize("/").last()
    name = probeset.take(probeset.lastIndexOf('.'))
    """
    wget ${address}
    """
}
