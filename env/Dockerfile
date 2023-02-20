FROM nfcore/base:2.1
LABEL \
    author = "Erik Fasterius" \
    description = "Temporary Dockerfile for nf-core/spatialtranscriptomics" \
    maintainer = "erik.fasterius@nbis.se"

# Copy and install the Conda environment
COPY environment.yml ./
RUN conda config --set channel_priority strict \
    && conda env update --name base --file environment.yml \
    && conda clean --all --force-pkgs-dirs --yes

# Installation of non-Conda packages
RUN Rscript -e 'devtools::install_github("MarcElosua/SPOTlight", ref="1f364a965ab275ac54a711e28c67554db25d547e")'
RUN Rscript -e 'remotes::install_github("JEFworks-Lab/STdeconvolve", ref="648d2303bb38002dc6e9a5c9f59487b44d8161d4")'

# Start Bash shell by default
CMD /bin/bash
