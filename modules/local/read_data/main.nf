//
// Read ST 10x Visium data with spatialdata_io and save to `SpatialData` file
//
process READ_DATA {

    tag "${meta.id}"
    label 'process_low'

    container "docker.io/erikfas/spatialvi"

    input:
    tuple val(meta), path("${meta.id}")
    val(hd_bin_size)

    output:
    tuple val(meta), path("sdata_raw.zarr"), emit: sdata_raw
    path("versions.yml")                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "The READ_DATA module does not support Conda/Mamba, please use Docker / Singularity / Podman instead."
    }

    """
    # Fix required directory structure

    # Prepare the --visium_hd flag and bin_size conditionally
    visiumHdFlag=""
    binSizeFlag=""
    if [ -d "${meta.id}/binned_outputs" ]; then
        visiumHdFlag="--visium_hd"
        binSizeFlag="--bin_size ${hd_bin_size}"
    fi

    # Set environment variables
    export XDG_CACHE_HOME="./.xdg_cache_home"
    export XDG_DATA_HOME="./.xdg_data_home"

    # Execute read data script
    read_data.py \\
        --SRCountDir "${meta.id}" \\
        --sampleID "${meta.id}" \\
        --output_sdata sdata_raw.zarr \\
        \$visiumHdFlag \\
        \$binSizeFlag

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spatialdata_io: \$(python -c "import spatialdata_io; print(spatialdata_io.__version__)")
    END_VERSIONS
    """
}
