//
// Download SpaceRanger probeset
//

process SPACERANGER_DOWNLOAD_PROBESET {

    tag "${name}"
    label "process_low"

    conda "conda-forge::gnu-wget=1.18"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--hed695b0_4' :
        'quay.io/biocontainers/gnu-wget:1.18--hed695b0_4' }"

    input:
    val(address)

    output:
    path("*.csv")       , emit: probeset
    path("versions.yml"), emit: versions

    script:
    probeset = address.tokenize("/").last()
    name = probeset.take(probeset.lastIndexOf('.'))
    """
    wget ${address}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gnu-wget: \$(wget --version | head -1 | cut -d ' ' -f 3)
    END_VERSIONS
    """
}
