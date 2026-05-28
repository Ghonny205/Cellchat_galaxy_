# 1. Arquitectura nativa
FROM condaforge/miniforge3:latest

# 2. Herramientas de sistema de Linux
RUN apt-get update && apt-get install -y \
    build-essential \
    libglpk-dev \
    gcc \
    g++ \
    gfortran \
    libblas-dev \
    liblapack-dev \
    pkg-config

# 3. Instalamos todo EXCEPTO NMF (que nos da la versión mala)
RUN mamba install -y -c conda-forge -c bioconda \
    r-base \
    r-seurat \
    r-optparse \
    r-devtools \
    r-remotes \
    r-igraph \
    r-circlize \
    r-ggpubr \
    r-ggnetwork \
    r-rstatix \
    r-nloptr \
    r-svglite \
    r-rspectra \
    bioconductor-complexheatmap \
    bioconductor-biocneighbors \
    bioconductor-biobase

# 4. Instalamos el NMF moderno y luego CellChat
RUN Rscript -e "install.packages('NMF', repos='http://cran.rstudio.com/'); \
                remotes::install_github('sqjin/CellChat', upgrade = 'never'); \
                if (!require('CellChat', character.only=TRUE)) quit(status=1)"