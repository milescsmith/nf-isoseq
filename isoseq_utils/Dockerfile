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
    conda env create --file /opt/environment.yml && \
    conda clean --all -y


CMD [ 'conda' ]