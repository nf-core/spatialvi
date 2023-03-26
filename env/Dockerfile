FROM nfcore/base:2.1
LABEL \
    author = "Christophe Avenel" \
    description = "Temporary Dockerfile for nf-core/spatialtranscriptomics" \
    maintainer = "christophe.avenel@it.uu.se"

RUN apt-get -y update; apt-get -y install -y --no-install-recommends \
    pandoc \
    pandoc-citeproc \
    curl \
    gdebi-core \
    && rm -rf /var/lib/apt/lists/*

# Install quarto
ARG QUARTO_VERSION="0.9.522"
RUN curl -o quarto-linux-amd64.deb -L https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb
RUN gdebi --non-interactive quarto-linux-amd64.deb

# Copy and install the Conda environment
COPY environment.yml ./
RUN conda install -c conda-forge mamba \
    && mamba env update --name base --file environment.yml \
    && mamba clean --all --force-pkgs-dirs --yes

# Start Bash shell by default
CMD /bin/bash
