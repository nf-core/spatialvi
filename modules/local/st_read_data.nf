//
// Read ST 10x visium and SC 10x data with Scanpy and save to `anndata` file
//
process ST_READ_DATA {

    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::scanpy=1.7.2 conda-forge::matplotlib=3.6.3 conda-forge::pandas=1.5.3"
    container "docker.io/cavenel/spatialtranscriptomics"

    input:
    tuple val (meta), path("${meta.id}/*")

    output:
    tuple val(meta), path("st_sdata.zarr"), emit: st_sdata_raw
    path("versions.yml")                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir "${meta.id}/spatial"
    mv  "${meta.id}/scalefactors_json.json" \\
        "${meta.id}/tissue_hires_image.png" \\
        "${meta.id}/tissue_lowres_image.png" \\
        "${meta.id}/tissue_positions.csv" \\
        "${meta.id}/spatial/"

    read_st_data.py \\
        --SRCountDir "${meta.id}" \\
        --output_sdata st_sdata.zarr

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """
}
