# nf-core/spatialtranscriptomics: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

Initial release of nf-core/spatialtranscriptomics, created with the
[nf-core](https://nf-co.re/) template. This marks the point at which the
pipeline development was moved to nf-core and NBIS. The pipeline has undergone
several iterations regarding its functionality and content; there are a
significant number of changes, of which not all are listed here. In summary, the
pipeline contains best-practice processing and analyses of pre- and post-Space
Ranger-processed data, including quality controls, normalisation, dimensionality
reduction, clustering, differential expression testing as well as output files
compatible with further downstream analyses and/or exploration in _e.g._
[TissUUmaps](https://tissuumaps.github.io/) or bespoke user code.

### `Added`

- Add MultiQC support for Space Ranger outputs [[#70](https://github.com/nf-core/spatialtranscriptomics/pull/70)]
- Use the QUARTONOTEBOOK nf-core module instead of local Quarto-based modules [[#68](https://github.com/nf-core/spatialtranscriptomics/pull/68)]
- Add a custom nf-core Quarto template for the downstream analysis reports [[#64](https://github.com/nf-core/spatialtranscriptomics/pull/64)]
- Allow input directories `fastq_dir` and `spaceranger_dir` to be specified as tar archives (`.tar.gz`)
- Add a check to make sure that there are spots left after filtering [[#46](https://github.com/nf-core/spatialtranscriptomics/issues/46)]
- Implement tests with nf-test [[#42](https://github.com/nf-core/spatialtranscriptomics/pull/42)]
- Replace custom code to download reference with `untar` module [[#44](https://github.com/nf-core/spatialtranscriptomics/pull/44)]
- Embed resources in quarto reports [[#43](https://github.com/nf-core/spatialtranscriptomics/pull/43)]
- Use a samplesheet for input specification [[#30](https://github.com/nf-core/spatialtranscriptomics/pull/30), [#31](https://github.com/nf-core/spatialtranscriptomics/pull/31) and [#45](https://github.com/nf-core/spatialtranscriptomics/pull/45)]
- Add Space Ranger pre-processing as an optional pipeline step using the `spaceranger` nf-core module [[#17](https://github.com/nf-core/spatialtranscriptomics/pull/17) and [#45](https://github.com/nf-core/spatialtranscriptomics/pull/45)]
- Add `env/` directory with pipeline-specific container and Conda environment specifications [[#17](https://github.com/nf-core/spatialtranscriptomics/pull/17) and [#28](https://github.com/nf-core/spatialtranscriptomics/pull/28)]
- Use a more standardised practice to find mitochondrial genes [[#30](https://github.com/nf-core/spatialtranscriptomics/pull/30)]
- Make pipeline output compatible with TissUUmaps [[#31](https://github.com/nf-core/spatialtranscriptomics/pull/31)]
- Add custom Quarto-based reports for all downstream processing [[#31](https://github.com/nf-core/spatialtranscriptomics/pull/31)]
- Embed resources in quarto reports [[#43](https://github.com/nf-core/spatialtranscriptomics/pull/43)]

### `Fixed`

- [#51](https://github.com/nf-core/spatialtranscriptomics/issues/51): Fix version export of `leidenalg` and `SpatialDE` Python modules
- [#38](https://github.com/nf-core/spatialtranscriptomics/issues/38): Specify manual alignment files in samplesheet
- [#20](https://github.com/nf-core/spatialtranscriptomics/issues/20) and [#22](https://github.com/nf-core/spatialtranscriptomics/issues/22): Add missing Groovy module
- [#53](https://github.com/nf-core/spatialtranscriptomics/pull/53): Use ensemble IDs as index in adata.var and fix related
  issue with SpatialDE

### `Dependencies`

Note, since the pipeline is using Nextflow DSL2, each process will be run
with its own [Biocontainer](https://biocontainers.pro/#/registry). This means
that on occasion it is entirely possible for the pipeline to be using different
versions of the same tool.

| Dependency  | Version |
| ----------- | ------- |
| `SpatialDE` | 1.1.3   |
| `leidenalg` | 0.9.1   |
| `python`    | 3.12.0  |
| `quarto`    | 1.3.302 |
| `scanpy`    | 1.9.3   |

### `Removed`

- Streamline pipeline for basic ST data processing; remove SC processing and deconvolution (for now) [[#31](https://github.com/nf-core/spatialtranscriptomics/pull/31)]

## v0.1.0 - 2023-03-31

Initial release of nf-core/spatialtranscriptomics, created with the
[nf-core](https://nf-co.re/) template by the Jackson Laboratory contributors
(see `README.md` for details).
