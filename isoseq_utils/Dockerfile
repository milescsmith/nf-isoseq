FROM ubuntu:bionic

LABEL author="Miles Smith <miles-smith@omrf.org>"
LABEL version="0.1.0"

COPY environment.yml /opt/environment.yml
ADD https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh /opt/miniconda.sh

ENV PATH /home/root/conda/bin:$PATH
RUN /bin/bash /opt/miniconda.sh -b -p /home/root/conda && \
    rm /opt/miniconda.sh && \
    /home/root/conda/bin/conda init && \
    conda update -n base -c defaults conda -y && \
    conda install -n base -c bioconda gmap pbccs==4.2.0 pbbam==1.3.0 lima==1.11.0 isoseq3==3.3.0 gmap samtools && \
    conda clean --all -y
