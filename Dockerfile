FROM vanallenlab/miniconda:3.6

WORKDIR /

RUN apt-get update \
    && apt-get -y install gcc build-essential gfortran git texlive libgomp1 libz-dev \
    && conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda install --yes htslib \
    && conda install --yes -c r r-openssl r-devtools r-ggplot2

ENV LD_LIBRARY_PATH /opt/conda/lib::$LD_LIBRARY_PATH

RUN Rscript -e "options(unzip = 'internal'); devtools::install_github('mskcc/facets@v0.5.14', build_vignettes = TRUE)"

RUN g++ -std=c++11 -I/opt/conda/pkgs/htslib-1.9-ha228f0b_7/include /opt/conda/lib/R/library/facets/extcode/snp-pileup.cpp -L/opt/conda/pkgs/htslib-1.9-ha228f0b_7/lib -lhts -Wl,-rpath=/opt/conda/pkgs/htslib-1.9-ha228f0b_7/lib -o snp-pileup

COPY facets.R Dockerfile facets.wdl README.md /
