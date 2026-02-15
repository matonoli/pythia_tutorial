FROM rootproject/root:6.28.00-ubuntu22.04

RUN apt-get update && apt-get install -y \
    git \
    python3-pip \
    wget \
    bzip2 \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install Miniforge (conda-forge distribution; avoids Anaconda pkgs/main & pkgs/r ToS)
RUN wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh && \
    bash /tmp/miniforge.sh -b -p /opt/conda && \
    rm /tmp/miniforge.sh && \
    ln -s /opt/conda/bin/conda /usr/local/bin/conda && \
    /opt/conda/bin/conda config --system --remove channels defaults || true && \
    /opt/conda/bin/conda config --system --add channels conda-forge && \
    /opt/conda/bin/conda config --system --set channel_priority strict

# Install Pythia and notebook deps from conda-forge only
# (Removed dask/distributed since they're not used in this repo's shared code.)
RUN /opt/conda/bin/conda install -y -c conda-forge \
    python=3.10.6 \
    pythia8 \
    jupyter \
    ipykernel \
    metakernel \
    && /opt/conda/bin/conda clean -afy

ENV PATH=/opt/conda/bin:$PATH
ENV LD_LIBRARY_PATH=/opt/conda/lib:$LD_LIBRARY_PATH
ENV PYTHIA8=/opt/conda

RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> /etc/bash.bashrc && \
    echo "conda activate base" >> /etc/bash.bashrc