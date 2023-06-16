//
// Read ST 10x visium and SC 10x data with Scanpy and save to `anndata` file
//
process ST_READ_DATA {

    tag "${meta.id}"
    label "process_low"

    conda "conda-forge::scanpy=1.7.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scanpy:1.7.2--pyhdfd78af_0' :
        'quay.io/biocontainers/scanpy:1.7.2--pyhdfd78af_0' }"

    input:
    tuple val (meta), path("${meta.id}")

    output:
    tuple val(meta), path("st_adata_raw.h5ad"), emit: st_raw
    path("versions.yml")                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    read_st_data.py \\
        --SRCountDir "${meta.id}" \\
        --outAnnData st_adata_raw.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """
}
