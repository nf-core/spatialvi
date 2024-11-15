//
// Merge per-sample SpatialData into a single SpatialData
//
process MERGE_SDATA {

    label 'process_low'
    container "docker.io/erikfas/spatialvi"

    input:
    path(sdata, stageAs: "?/*")

    output:
    path("merged_sdata.zarr"), emit: sdata
    path("versions.yml")     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "The MERGE_SDATA module does not support Conda/Mamba, please use Docker / Singularity / Podman instead."
    }
    """
    # Set environment variables
    export XDG_CACHE_HOME="./.xdg_cache_home"
    export XDG_DATA_HOME="./.xdg_data_home"

    # Execute script
    merge_sdata.py \\
        ${sdata} \\
        merged_sdata.zarr

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spatialdata_io: \$(python -c "import spatialdata_io; print(spatialdata_io.__version__)")
    END_VERSIONS
    """
}
