# First stage: multi-platform Quarto image
FROM jdutant/quarto-minimal:1.3.313 as quarto

# Second stage: multi-platform Mamba image
FROM condaforge/mambaforge:23.11.0-0

# Copy Quarto installation from first stage and add to PATH
COPY --from=quarto /opt/quarto /opt/quarto
ENV PATH="${PATH}:/opt/quarto/bin"

# Install packages using Mamba; also remove static libraries, python bytecode
# files and javascript source maps that are not required for execution
COPY environment.yml ./
RUN mamba env update --name base --file environment.yml \
    && mamba clean --all --force-pkgs-dirs --yes \
    && find /opt/conda -follow -type f -name '*.a' -delete \
    && find /opt/conda -follow -type f -name '*.pyc' -delete \
    && find /opt/conda -follow -type f -name '*.js.map' -delete

CMD /bin/bash

LABEL \
    authors = "Erik Fasterius, Christophe Avenel" \
    description = "Dockerfile for nf-core/spatialtranscriptomics report modules"
