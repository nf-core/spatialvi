# ![nf-core/spatialtranscriptomics](docs/images/nf-core-spatialtranscriptomics_logo_light.png#gh-light-mode-only) ![nf-core/spatialtranscriptomics](docs/images/nf-core-spatialtranscriptomics_logo_dark.png#gh-dark-mode-only)

[![GitHub Actions CI Status](https://github.com/nf-core/spatialtranscriptomics/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/spatialtranscriptomics/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/nf-core/spatialtranscriptomics/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/spatialtranscriptomics/actions?query=workflow%3A%22nf-core+linting%22)
[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/st/results)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23st-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/st)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)
[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/spatialtranscriptomics** is a bioinformatics best-practice analysis pipeline for Spatial Transcriptomics Integrated Analysis.
The pipeline for processing spatially-resolved gene counts with spatial coordinates, image data, and optionally single cell RNA-seq data, designed for 10x genomics visium and single cell transcriptomics. Specifically, input data can be output of 10x SpaceRanger and CellRanger.

The are numerous methods for ST data analysis, and this research area is rapidly developing. The pipeline may be useful to a community working in the area of ST. The pipeline that performs a set of analyses, not limited to but including quality control, normalization, deconvolution of spots into cell types and topics, resolution enhancement, clustering, selection, spatially-variable genes, etc. The software is made of R and Python with numerous packages containerized.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

<!-- TODO nf-core: Add full-sized test dataset and amend the paragraph below if applicable -->
On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/st/results).

## Pipeline summary

1. Normalization and Quality Control ([`scran`](https://doi.org/doi:10.18129/B9.bioc.scran), [`scanpy`](https://github.com/theislab/scanpy))
2. Spots cell types and cell topics deconvolution ([`STdeconvolve`](https://jef.works/STdeconvolve/), [`SPOTlight`](https://github.com/MarcElosua/SPOTlight))
3. Spot resolution enhancement ([`BayesSpace`](https://github.com/edward130603/BayesSpace))
4. Dimensionality reduction and projection ([`scanpy`](https://github.com/theislab/scanpy), [`Seurat`](https://satijalab.org/seurat/))
5. Integration with scRNA-seq data ([`scanorama`](https://github.com/brianhie/scanorama), [`BBKNN`](https://github.com/Teichlab/bbknn)) 
6. Label transfer from scRNA-seq data ([`Seurat`](https://satijalab.org/seurat/))
7. Clustering of the spots ([`scanpy Leiden`](https://arxiv.org/abs/1810.08473), [`BayesSpace`](https://github.com/edward130603/BayesSpace))
8. Visualization of clusters and features in spatial coordinates and 2D projection layout ([`scanpy`](https://github.com/theislab/scanpy), [`Seurat`](https://satijalab.org/seurat/))
9. Identification of spatially variable features ([`SpatialDE`](https://github.com/Teichlab/SpatialDE))


The pipeline combines multiple tools, toolkits and platforms:

+ [`Bioconductor`](https://www.bioconductor.org/) - software resource for the analysis of genomic data. Based primarily on the statistical R programming language.
+ [`Seurat`](https://satijalab.org/seurat/) - R toolkit for single cell genomics.
+ [`scran`](https://doi.org/doi:10.18129/B9.bioc.scran) - R package implements miscellaneous functions for analysis and interpretation of single-cell RNA-seq data.
+ [`SpatialExperiment`](https://doi.org/doi:10.18129/B9.bioc.SpatialExperiment) - R package for storing, retrieving spatial coordinates and gene expression.
+ [`Giotto`](https://rubd.github.io/Giotto_site/) - (R package) a toolbox for integrative analysis and visualization of spatial expression data.
+ [`reticulate`](https://github.com/rstudio/reticulate/) - comprehensive set of tools for interoperability between Python and R.
+ [`STdeconvolve`](https://jef.works/STdeconvolve/) - R implementation of LDA-based cell-topics deconvolution of spots.
+ [`SPOTlight`](https://github.com/MarcElosua/SPOTlight) - R implementation of NMF-based cell-types deconvolution of spots.
+ [`BayesSpace`](https://github.com/edward130603/BayesSpace) - R package for spatially-aware clustering and resolution enhancement.
+ [`SpatialDE`](https://github.com/Teichlab/SpatialDE) - Python package for identification of spatially variable genes.
+ [`scanorama`](https://github.com/brianhie/scanorama) - Python package that enables batch-correction and integration of heterogeneous scRNA-seq datasets.
+ [`BBKNN`](https://github.com/Teichlab/bbknn) - (Python package) batch effect removal tool.
+ [`scanpy`](https://github.com/theislab/scanpy) - scalable Python toolkit for analyzing single-cell gene expression data.
+ [`anndata`](https://github.com/theislab/anndata) - a Python package for handling annotated data matrices in memory and on disk.


## Quick Start

> Note: As a temporary measure the singularity containers necessary to run this pipeline were uploaded to https://doi.org/10.5281/zenodo.6266244. Download the two containers and edit the `nextflow.config` parameters `container_python` and `container_r`.

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```console
    nextflow run nf-core/spatialtranscriptomics -profile test,YOURPROFILE
    ```

    Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

    > * The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
    > * Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
    > * If you are using `singularity` and are persistently observing issues downloading Singularity images directly due to timeout or network issues, then you can use the `--singularity_pull_docker_container` parameter to pull and convert the Docker image instead. Alternatively, you can use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
    > * If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

    <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

    ```console
    nextflow run nf-core/spatialtranscriptomics -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> --input samplesheet.csv --genome GRCh37
    ```

## Documentation

The nf-core/spatialtranscriptomics pipeline comes with documentation about the pipeline [usage](https://nf-co.re/st/usage), [parameters](https://nf-co.re/st/parameters) and [output](https://nf-co.re/st/output).

## Credits

nf-core/spatialtranscriptomics was originally developed by The Jackson Laboratory as a part of [NCI PDXNet PDCCC](https://www.pdxnetwork.org/pdccc) (PDX Data Commons and Coordination Center). Original authors:

+ [Sergii Domanskyi](https://github.com/sdomanskyi)
+ [Jeffrey Chuang]()
+ [Anuj Srivastava]()

The pipeline is being further developed in collaboration with the [National Genomics Infastructure](https://ngisweden.scilifelab.se/) within [SciLifeLab](https://scilifelab.se/).

<!-- We thank the following people for their extensive assistance in the development of this pipeline: -->

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#st` channel](https://nfcore.slack.com/channels/st) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/spatialtranscriptomics for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->
An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
