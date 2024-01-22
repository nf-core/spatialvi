//
// Quality controls and filtering
//
process ST_QUALITY_CONTROLS {

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
    tuple val(meta), path(st_adata_raw)

    output:
    tuple val(meta), path("st_adata_filtered.h5ad")  , emit: st_adata_filtered
    tuple val(meta), path("st_quality_controls.html"), emit: html
    path("versions.yml")                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    quarto render ${report} \
        -P input_adata_raw:${st_adata_raw} \
        -P min_counts:${params.st_preprocess_min_counts} \
        -P min_genes:${params.st_preprocess_min_genes} \
        -P min_spots:${params.st_preprocess_min_spots} \
        -P output_adata_filtered:st_adata_filtered.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto -v)
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """
}
