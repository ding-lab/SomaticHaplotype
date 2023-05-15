FROM rocker/tidyverse:4.1.0

RUN apt-get update && apt-get install -y --no-install-recommends \
    libglu1

# Install some Bioconductor packages dependencies
RUN Rscript -e "options(warn = 2); BiocManager::install('IRanges', update = FALSE)"
RUN Rscript -e "options(warn = 2); BiocManager::install('limma', update = FALSE)"

# R packages
RUN install2.r --error --deps TRUE NORMT3

# ref = 341eb77105e7efd2654b4f112578648584936e06 is latest greenelab/TDM commit (retrieved 2021-05-28)
RUN Rscript -e "options(warn = 2); remotes::install_github( \
    'genome/bmm')"

# ref = 08ed6b54e4efe5249107cb335cd8e169657cbc44 is wgmao/PLIER commit corresponding to v0.1.6 (retrieved 2021-11-09)
RUN Rscript -e "options(warn = 2); remotes::install_github( \
    'genome/sciClone')"
