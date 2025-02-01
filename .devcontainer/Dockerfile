FROM rootproject/root:6.28.00-ubuntu22.04

# # Run the following commands as super user (root):
# USER root

RUN apt-get update && apt-get install -y git python3-pip

# # Install Conda

RUN apt-get install -y wget bzip2 && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh && \
    /opt/conda/bin/conda init && \
    ln -s /opt/conda/bin/conda /usr/bin/conda

# # Install Pythia and other dependencies via Conda
RUN /opt/conda/bin/conda install -c conda-forge -y python=3.10.6 \
    pythia8 \
    jupyter \
    ipykernel \
    metakernel \
    dask \
    distributed && \
    /opt/conda/bin/conda clean -afy