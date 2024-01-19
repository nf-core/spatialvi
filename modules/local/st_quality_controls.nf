//
// Spatial data pre-processing
//
process ST_QUALITY_CONTROLS {

    // TODO: Add a better description
    // TODO: Update Conda directive when Quarto/Pandoc works on ARM64

    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::quarto=1.3.353 conda-forge::scanpy=1.9.3 conda-forge::papermill=2.3.4 conda-forge::jupyter=1.0.0"
    container "docker.io/erikfas/spatialtranscriptomics"

    // Exit if running this module with -profile conda / -profile mamba on ARM64
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        architecture = System.getProperty("os.arch")
        if (architecture == "arm64" || architecture == "aarch64") {
            exit 1, "The ST_QUALITY_CONTROLS module does not support Conda on ARM64. Please use Docker / Singularity / Podman instead."
        }
    }

    input:
    path(report)
    path(report_template)
    tuple val(meta), path(st_raw, stageAs: "adata_raw.h5ad")

    output:
    tuple val(meta), path("st_adata_norm.h5ad")           , emit: st_data_norm
    tuple val(meta), path("st_quality_controls.html") , emit: html
    path("versions.yml")                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    quarto render ${report} \
        --output st_quality_controls.html \
        -P rawAdata:${st_raw} \
        -P minCounts:${params.st_preprocess_min_counts} \
        -P minGenes:${params.st_preprocess_min_genes} \
        -P minCells:${params.st_preprocess_min_spots} \
        -P nameDataNorm:st_adata_norm.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto -v)
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """
}
