//
// Spatial data pre-processing
//
process ST_QC_AND_NORMALISATION {

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
            exit 1, "The ST_QC_AND_NORMALISATION module does not support Conda on ARM64. Please use Docker / Singularity / Podman instead."
        }
    }

    input:
    path(report)
    tuple val(meta), path(st_raw, stageAs: "adata_raw.h5ad")

    output:
    tuple val(meta), path("st_adata_norm.h5ad")           , emit: st_data_norm
    tuple val(meta), path("st_adata_plain.h5ad")          , emit: st_data_plain
    tuple val(meta), path("st_qc_and_normalisation.html") , emit: html
    path("versions.yml")                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    quarto render ${report} \
        --output st_qc_and_normalisation.html \
        -P rawAdata:${st_raw} \
        -P pltFigSize:${params.st_preprocess_fig_size} \
        -P minCounts:${params.st_preprocess_min_counts} \
        -P minGenes:${params.st_preprocess_min_genes} \
        -P minCells:${params.st_preprocess_min_cells} \
        -P histplotQCmaxTotalCounts:${params.st_preprocess_hist_qc_max_total_counts} \
        -P histplotQCminGeneCounts:${params.st_preprocess_hist_qc_min_gene_counts} \
        -P histplotQCbins:${params.st_preprocess_hist_qc_bins} \
        -P nameDataPlain:st_adata_plain.h5ad \
        -P nameDataNorm:st_adata_norm.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quarto: \$(quarto -v)
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """
}
