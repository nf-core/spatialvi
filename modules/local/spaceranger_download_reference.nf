//
// Download genome reference
//

process SPACERANGER_DOWNLOAD_REFERENCE {

    tag "${name}"
    label "process_low"

    conda "conda-forge::gnu-wget=1.18"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--hed695b0_4' :
        'quay.io/biocontainers/gnu-wget:1.18--hed695b0_4' }"

    input:
    val(address)

    output:
    path "${name}", type: "dir", emit: reference
    path "versions.yml"        , emit: versions

    script:
    reference = address.tokenize("/").last()
    name = reference.split("\\.", 2)[0]
    """
    wget ${address}
    tar -xzvf ${reference}
    rm ${reference}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gnu-wget: \$(wget --version | head -1 | cut -d ' ' -f 3)
    END_VERSIONS
    """
}
