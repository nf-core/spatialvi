<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-spatialtranscriptomics_logo_dark.png">
    <img alt="nf-core/spatialtranscriptomics" src="docs/images/nf-core-spatialtranscriptomics_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/spatialtranscriptomics/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/spatialtranscriptomics/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/spatialtranscriptomics/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/spatialtranscriptomics/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/spatialtranscriptomics/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/spatialtranscriptomics)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23spatialtranscriptomics-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/spatialtranscriptomics)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/spatialtranscriptomics** is a bioinformatics analysis pipeline for
Spatial Transcriptomics. It can process and analyse 10X spatial data either
directly from raw data by running [Space Ranger](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/what-is-space-ranger)
or data already processed by Space Ranger. The pipeline currently consists of the
following steps:

0. Raw data processing with Space Ranger (optional)
1. Quality controls and filtering
2. Normalisation
3. Dimensionality reduction and clustering
4. Differential gene expression testing

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool
to run tasks across multiple compute infrastructures in a very portable manner.
It uses Docker/Singularity containers making installation trivial and results
highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html)
implementation of this pipeline uses one container per process which makes it
much easier to maintain and update software dependencies. Where possible, these
processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules)
in order to make them available to all nf-core pipelines, and to everyone within
the Nextflow community!

On release, automated continuous integration tests run the pipeline on a
full-sized dataset on the AWS cloud infrastructure. This ensures that the
pipeline runs on AWS, has sensible resource allocation defaults set to run on
real-world datasets, and permits the persistent storage of results to benchmark
between pipeline releases and other analysis sources. The results obtained from
the full-sized test can be viewed on the [nf-core website](https://nf-co.re/spatialtranscriptomics/results).

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

You can run the pipeline using:

```bash
nextflow run nf-core/spatialtranscriptomics \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/spatialtranscriptomics/usage) and the [parameter documentation](https://nf-co.re/spatialtranscriptomics/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/spatialtranscriptomics/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/spatialtranscriptomics/output).

## Credits

nf-core/spatialtranscriptomics was originally developed by the Jackson
Laboratory<sup>1</sup>, up to the [0.1.0](https://github.com/nf-core/spatialtranscriptomics/releases/tag/0.1.0)
tag. It was further developed in a collaboration between the [National
Bioinformatics Infrastructure Sweden](https://nbis.se/) and [National Genomics
Infrastructure](https://ngisweden.scilifelab.se/) within [SciLifeLab](https://scilifelab.se/);
it is currently developed and maintained by [Erik Fasterius](https://github.com/fasterius)
and [Christophe Avenel](https://github.com/cavenel).

Many thanks to others who have helped out along the way too, especially [Gregor
Sturm](https://github.com/grst)!

_<sup>1</sup> Supported by grants from the US National Institutes of Health
[U24CA224067](https://reporter.nih.gov/project-details/10261367) and
[U54AG075941](https://reporter.nih.gov/project-details/10376627). Original
authors [Dr. Sergii Domanskyi](https://github.com/sdomanskyi), Prof. Jeffrey
Chuang and Dr. Anuj Srivastava._

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#spatialtranscriptomics` channel](https://nfcore.slack.com/channels/spatialtranscriptomics) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/spatialtranscriptomics for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
