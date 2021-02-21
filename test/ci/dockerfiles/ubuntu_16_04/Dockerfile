FROM ubuntu:16.04
MAINTAINER David S. Kammer <dkammer@ethz.ch>

ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8

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
    python3-sphinx \
    python3-sphinx-rtd-theme \
    && rm -rf /var/lib/apt/lists/*
