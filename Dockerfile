FROM vanallenlab/miniconda:3.6

WORKDIR /

RUN apt-get update \
    && apt-get -y install gcc git \
    && conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda install --yes htslib gcc \
    && git clone --branch 'v0.5.14' https://github.com/mskcc/facets.git

ENV LD_LIBRARY_PATH /opt/conda/pkgs/libdeflate-1.0-h470a237_0/lib::$LD_LIBRARY_PATH

RUN g++ -std=c++11 -I/opt/conda/pkgs/htslib-1.9-hc238db4_4/include /facets/inst/extcode/snp-pileup.cpp -L/opt/conda/pkgs/htslib-1.9-hc238db4_4/lib -lhts -Wl,-rpath=/opt/conda/pkgs/htslib-1.9-hc238db4_4/lib -o snp-pileup