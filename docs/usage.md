# nf-core/spatialtranscriptomics: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/spatialtranscriptomics/usage](https://nf-co.re/spatialtranscriptomics/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Samplesheet input

You will need to create a samplesheet with information about the samples you
would like to analyse before running the pipeline. It has to be a comma-separated file as described
in the examples below and depends on the input data type. Use this parameter to specify its location.

```bash
--input '[path to samplesheet file]'
```

There are two types of samplesheets that the pipeline can handle: those
specifying _raw data_ (to be analysed by Space Ranger) and _processed data_
(_i.e._ already analysed by Space Ranger). The workflow will automatically
detect the samplesheet type and run the appropriate analysis steps. The two
types of samplesheet are described in the following sections.

### Raw spatial data

This section describes samplesheets for processing _raw spatial data_ yet to be analysed with Space Ranger.

Here is an example of a typical samplesheet for analysing FFPE or fresh frozen (FF) data with bright field microscopy
imagery:

```no-highlight
sample,fastq_dir,image,slide,area
SAMPLE_1,fastqs_1/,hires_1.png,V11J26,B1
SAMPLE_2,fastqs_2/,hires_2.png,V11J26,B1
```

You may also supply a compressed tarball containing the FASTQ files in lieu of a
directory path:

```no-highlight
sample,fastq_dir,image,slide,area
SAMPLE_1,fastqs_1.tar.gz,hires_1.png,V11J26,B1
SAMPLE_2,fastqs_2.tar.gz,hires_2.png,V11J26,B1
```

For Cytassist samples, the `image` column gets replaced with the `cytaimage` column:

```no-highlight
sample,fastq_dir,cytaimage,slide,area
SAMPLE_1,fastqs_1/,cytassist_1.tif,V11J26,B1
SAMPLE_2,fastqs_2/,cytassist_2.tif,V11J26,B1
```

Depending on the experimental setup, (additional) colour composite fluorescence images or dark background
fluorescence images can be supplied using the `colorizedimage` or `darkimage` columns, respectively.

Please refer to the following table for an overview of all supported columns:

| Column             | Description                                                                                                           |
| ------------------ | --------------------------------------------------------------------------------------------------------------------- |
| `sample`           | Unique sample identifier. MUST match the prefix of the fastq files                                                    |
| `fastq_dir`        | Path to directory where the sample FASTQ files are stored. May be a `.tar.gz` file instead of a directory.            |
| `image`            | Brightfield microscopy image                                                                                          |
| `cytaimage`        | Brightfield tissue image captured with Cytassist device                                                               |
| `colorizedimage`   | A colour composite of one or more fluorescence image channels saved as a single-page, single-file colour TIFF or JPEG |
| `darkimage`        | Dark background fluorescence microscopy image                                                                         |
| `slide`            | The Visium slide ID used for the sequencing.                                                                          |
| `area`             | Which slide area contains the tissue sample.                                                                          |
| `manual_alignment` | Path to the manual alignment file (optional)                                                                          |
| `slidefile`        | Slide specification as JSON. Overrides `slide` and `area` if specified. (optional)                                    |

> [!NOTE]
> - You need to specify _at least one_ of `image`, `cytaimage`, `darkimage`,
>   `colorizedimage`. Most commonly, you'll specify `image` for bright field
>   microscopy data, or `cytaimage` for tissue scans generated with the 10x
>   Cyatassist device. Please refer to the [Space Ranger documentation](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/what-is-space-ranger),
>   how multiple image types can be combined.
> - The `manual_alignment` column is only required for samples for which a
>   manual alignment file is needed and can be ignored if you're using automatic
>   alignment.

If you are unsure, please see the Visium documentation for details regarding the
different variants of [FASTQ directory structures](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/using/fastq-input)
and [slide parameters](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/using/slide-info)
appropriate for your samples.

### Processed data

If your data has already been processed by Space Ranger and you are only
interested in running downstream steps, the samplesheet looks as follows:

```no-highlight
sample,spaceranger_dir
SAMPLE_1,results/SAMPLE_1/outs
SAMPLE_2,results/SAMPLE_2/outs
```

You may alternatively supply a compressed tarball containing the Space Ranger output:

```no-highlight
sample,spaceranger_dir
SAMPLE_1,outs.tar.gz
SAMPLE_2,outs.tar.gz
```

| Column            | Description                                                                               |
| ----------------- | ----------------------------------------------------------------------------------------- |
| `sample`          | Unique sample identifier.                                                                 |
| `spaceranger_dir` | Output directory generated by spaceranger. May be a `.tar.gz` file instead of a directory |

The Space Ranger output directory is typically called `outs` and contains both
gene expression matrices as well as spatial information.

## Space Ranger

The pipeline exposes several of Space Ranger's parameters when executing with
raw spatial data. Space Ranger requieres a lot of memory
(64 GB) and several threads (8) to be able to run. You can find the Space Ranger
documentation at the [10X website](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/what-is-space-ranger).

You are only able to run Space Ranger on the [officially supported organisms](https://support.10xgenomics.com/spatial-gene-expression/software/downloads/latest):
human and mouse. If you have already downloaded a reference you may supply the
path to its directory (or another link from the 10X website above) using the
`--spaceranger_reference` parameter, otherwise the pipeline will download the
default human reference for you automatically.

> [!NOTE]
> For FFPE and Cytassist experiments, you need to manually supply the appropriate probset using the `--spaceranger_probeset` parameter
> Please refer to the [Spaceranger Downloads page](https://support.10xgenomics.com/spatial-gene-expression/software/downloads/latest)
> to obtain the correct probeset.

## Analysis options

The pipeline uses Python and the `scverse` tools to do the downstream analysis
(quality control, filtering, clustering, spatial differential equations).

### Parameters for Quality Control and Filtering:

The following parameters are exposed for preprocessing:

- `--st_preprocess_min_counts`: Minimum number of counts for a spot to be considered in the analysis.
- `--st_preprocess_min_genes`: Minimum number of genes expressed in a spot for the spot to be considered.
- `--st_preprocess_min_cells`: Minimum number of spots expressing a gene for the gene to be considered.
- `--st_preprocess_fig_size`: The figure size for the plots generated during preprocessing (_e.g._, quality control plots).
- `--st_preprocess_hist_qc_max_total_counts`: Maximum total counts for the histogram plot in quality control.
- `--st_preprocess_hist_qc_min_gene_counts`: Minimum gene counts for the histogram plot in quality control.
- `--st_preprocess_hist_qc_bins`: Number of bins for the histogram plot in quality control.

### Parameters for Clustering :

- `--st_cluster_resolution`: Resolution parameter for the clustering algorithm, controlling granularity.

### Parameters for Spatial Differential Expression :

- `st_spatial_de_ncols`: Number of columns in the output figure.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run \
    nf-core/spatialtranscriptomics \
    --input <SAMPLESHEET> \
    --outdir <OUTDIR> \
    -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/spatialtranscriptomics -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: '<SAMPLESHEET>'
outdir: '<OUTDIR>'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/spatialtranscriptomics
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/spatialtranscriptomics releases page](https://github.com/nf-core/spatialtranscriptomics/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

To further assist in reproducibility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

## Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!INFO]
> We highly recommend the use of Docker or Singularity containers for full
> pipeline reproducibility, however when this is not possible, Conda is
> partially supported. Please note that Conda is not at all supported for Space
> Ranger processing, and only supported on non-ARM64 architectures for analyses
> downstream of Space Ranger.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
