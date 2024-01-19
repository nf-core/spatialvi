//
// Clustering etc.
//
process ST_CLUSTERING {

    // TODO: Add a better description
    // TODO: Update Conda directive when Quarto/Pandoc works on ARM64

    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::quarto=1.3.353 conda-forge::scanpy=1.9.3 conda-forge::papermill=2.3.4 conda-forge::jupyter=1.0.0 conda-forge::leidenalg=0.9.1"
    container "docker.io/erikfas/spatialtranscriptomics"

    // Exit if running this module with -profile conda / -profile mamba on ARM64
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        architecture = System.getProperty("os.arch")
        if (architecture == "arm64" || architecture == "aarch64") {
            exit 1, "The ST_CLUSTERING module does not support Conda on ARM64. Please use Docker / Singularity / Podman instead."
        }
    }

    input:
    path(report)
    path(report_template)
    tuple val(meta), path(st_adata_norm, stageAs: "adata_norm.h5ad")

    output:
    tuple val(meta), path("st_adata_processed.h5ad"), emit: st_adata_processed
    tuple val(meta), path("st_clustering.html")     , emit: html
    path("versions.yml")                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    quarto render ${report} \
        --output "st_clustering.html" \
        -P fileNameST:${st_adata_norm} \
        -P resolution:${params.st_cluster_resolution} \
        -P nHighlyVariableGenes:${params.st_preprocess_n_hvgs} \
        -P saveFileST:st_adata_processed.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto -v)
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """
}
