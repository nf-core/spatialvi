# nf-core/spatialvi: Output

## Introduction

This document describes the output produced by the pipeline. Most of the output
is contained within HTML reports created with [Quarto](https://quarto.org/), but
there are also other files which you can either take and analyse further by
yourself or explore interactively with _e.g._ [TissUUmaps](https://tissuumaps.github.io/).

The directories listed below will be created in the results directory after the
pipeline has finished. Results for individual samples will be created in
subdirectories following the `<OUTDIR>/<SAMPLE>/` structure. All paths are
relative to the top-level results directory.

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes
data using the following steps:

- [Space Ranger](#space-ranger)
- [FastQC](#fastqc)
- [Per-sample data](#per-sample-data)
- [Per-sample reports](#per-sample-reports)
- [Integration](#integration)
- [Workflow reporting](#workflow-reporting)
  - [Pipeline information](#pipeline-information) - Report metrics generated
    during the workflow execution

## Space Ranger

<details markdown="1">
<summary>Output files</summary>

- `<SAMPLE>/spaceranger/`
  - `outs/spatial/tissue_[hi/low]res_image.png`: High and low resolution images.
  - `outs/spatial/tissue_positions.csv`: Spot barcodes and their array
    positions.
  - `outs/spatial/scalefactors_json.json`: Scale conversion factors for the
    spots.
  - `outs/filtered_feature_bc_matrix/barcodes.tsv.gz`: List of barcode IDs.
  - `outs/filtered_feature_bc_matrix/features.tsv.gz`: List of feature IDs.
  - `outs/filtered_feature_bc_matrix/matrix.mtx.gz`: Matrix of UMIs, barcodes
    and features.
  - `outs/web_summary.html`: Interactive summary report from Space Ranger.
  - `outs/cloupe.cloupe`: File for visualization in 10X Loupe Browser.
</details>

All files produced by Space Ranger are currently published as output of this
pipeline, regardless if they're being used downstream or not; you can find more
information about these files at the [10X website](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/overview).

## FastQC

<details markdown="1">
<summary>Output files</summary>

- `<SAMPLE>/fastqc/`
  - `*_fastqc.html`: FastQC report.
  - `*_fastqc.zip`: FastQC data archive.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives
general quality metrics about your sequenced reads. It provides information
about the quality score distribution across your reads and per base sequence
content among other metrics.

## Per-sample data

<details markdown="1">
<summary>Output files</summary>

- `<SAMPLE>/data/`
  - `<SAMPLE>-raw.h5ad`: Raw (unprocessed) data in AnnData format, as extracted
    from the SpatialData object before any downstream analysis.
  - `<SAMPLE>.h5ad`: Processed data in AnnData format, including QC metrics,
    clustering, differential expression, and spatial analysis results.
  - `<SAMPLE>_svg.csv`: List of spatially variable genes.

</details>

Data in `.h5ad` formats as processed by the pipeline, which can be
used for further downstream analyses if desired; unprocessed data is also
present in these files. It can also be used by the [TissUUmaps](https://tissuumaps.github.io/)
browser-based tool for visualisation and exploration, allowing you to delve into
the data in an interactive way. The list of spatially variable genes are added
as a convenience if you want to explore them in _e.g._ Excel.

## Per-sample reports

<details markdown="1">
<summary>Output files</summary>

- `<SAMPLE>/reports/`
  - `report-<SAMPLE>.html`: Analysis report including quality controls, clustering,
    differential expression, and spatial analysis results.
  - `report.qmd`: Quarto source notebook of the report.
  - `report.yml`: Parameters used by the report.
  - `_extensions/`: Quarto nf-core extension, common to all reports.
- `<SAMPLE>/plots/`
  - `*.png`: Individual plots from the per-sample report, saved as PNG files for
    easy reuse.

</details>

The per-sample reports contain quality control metrics, clustering results,
differential expression analysis, neighbourhood enrichment, cluster interaction
matrices, and spatially variable gene analysis.

## Integration

Integration outputs are only produced when `--integrate_data` is enabled (default)
and `--skip_downstream` is not set. This merges all per-sample data and applies
batch correction using the selected integration method (Harmony or Scanorama).

<details markdown="1">
<summary>Output files</summary>

- `integration/data/`
  - `merged.zarr`: Merged SpatialData object containing all samples (before
    integration).
  - `<METHOD>.zarr`: Integrated SpatialData object with batch-corrected results.
  - `<METHOD>.h5ad`: Integrated AnnData object with batch-corrected embeddings
    and clustering.
- `integration/reports/`
  - `report-integrated.html`: Integration report with cross-sample comparisons,
    batch mixing assessment, and integrated clustering results.
  - `report-integrated.qmd`: Quarto source notebook of the integration report.
  - `report-integrated.yml`: Parameters used by the integration report.
  - `_extensions/`: Quarto nf-core extension.
- `integration/plots/`
  - `*.png`: Individual plots from the integration report.

</details>

## Workflow reporting

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.
- `multiqc/`
  - Report generated by MultiQC: `multiqc_report.html`. Includes FastQC results,
    Space Ranger metrics, per-sample filtering statistics (spots and genes
    removed at each QC step), and software versions.
  - Data and plots generated by MultiQC: `multiqc_data/` and `multiqc_plots/`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent
functionality for generating various reports relevant to the running and
execution of the pipeline. This will allow you to troubleshoot errors with the
running of the pipeline, and also provide you with other information such as
launch commands, run times and resource usage.
