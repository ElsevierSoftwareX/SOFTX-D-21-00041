FROM ubuntu:18.04
MAINTAINER David S. Kammer <dkammer@ethz.ch>

ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8

## for apt to be noninteractive
ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true

RUN apt-get -qq update && apt-get -y -qq install \
    cmake \
    g++ \
    git \
    libgsl-dev \
    libatlas-base-dev \
    libfftw3-dev \
    libopenmpi-dev \
    openmpi-bin \
    xsltproc \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

COPY .openmpi /root/.openmpi

RUN python3 -m pip install sphinx sphinx_rtd_theme