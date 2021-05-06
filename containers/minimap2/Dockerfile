FROM ubuntu:hirsute-20210422

ARG MINIMAP2_VERSION="2.18"

LABEL base.image="ubuntu:hirsuite"
LABEL software="Minimap2"
LABEL software.version=${MINIMAP2_VERSION}}
LABEL website="https://github.com/lh3/minimap2"
LABEL maintainer1="Miles Smith"
LABEL maintainer1.email="miles-smith@omrf.org"

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update --fix-missing && \
    apt-get install --no-install-recommends -y \
        bzip2 && \
    apt-get clean && \
    rm -rf /tmp/downloaded_packages/* && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /opt
ADD https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2 /opt/minimap2.tar.bz2
RUN tar xfj /opt/minimap2.tar.bz2
ENV PATH="${PATH}:/opt/minimap2-${MINIMAP2_VERSION}_x64-linux"

ENTRYPOINT [ "minimap2" ]