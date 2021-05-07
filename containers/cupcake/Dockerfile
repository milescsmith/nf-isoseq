FROM python:3.8-slim as builder

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update --fix-missing && \
    apt-get install --no-install-recommends -y \
        gcc \
        git \
        python3-dev && \ 
    apt-get clean && \
    rm -rf /tmp/downloaded_packages/* && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /opt
RUN pip wheel --wheel-dir /opt --use-pep517 git+https://github.com/milescsmith/cDNA_cupcake#egg=cupcake

FROM python:3.8-slim as runner
WORKDIR /opt
COPY --from=builder /opt/cupcake-21.2.0-cp38-cp38-manylinux_2_28_x86_64.whl .
RUN pip install --no-cache-dir /opt/cupcake-21.2.0-cp38-cp38-manylinux_2_28_x86_64.whl && \
  rm /opt/cupcake-21.2.0-cp38-cp38-manylinux_2_28_x86_64.whl