FROM vanallenlab/miniconda:3.6

WORKDIR /

COPY requirements.txt /

RUN conda install --yes -c conda-forge --file requirements.txt

COPY Dockerfile /
COPY README.md /
COPY infer_wgd.py /
