//
// Read ST 10x visium and SC 10x data with spatialdata_io and save to `SpatialData` file
// TODO: Add Visium HD support
//
process READ_DATA {

    tag "${meta.id}"
    label 'process_low'

    container "docker.io/erikfas/spatialvi"

    input:
    tuple val (meta), path("${meta.id}")

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

    # Prepare the --visium_hd flag conditionally
    visiumHdFlag=""
    if [ -d "${meta.id}/binned_outputs" ]; then
        visiumHdFlag="--visium_hd"
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
        --bin_size 8

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spatialdata_io: \$(python -c "import spatialdata_io; print(spatialdata_io.__version__)")
    END_VERSIONS
    """
}
