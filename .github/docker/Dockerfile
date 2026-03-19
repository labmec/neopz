FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install necessary packages for building and running the application
RUN apt-get update && apt-get upgrade -y && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    g++ \
    gfortran \
    ccache \
    cmake \
    ninja-build \
    git \
    wget \
    ca-certificates \
    gnupg \
    liblapack-dev \
    libblas-dev \
    libboost-all-dev \
    python3 && \
    wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2023.PUB -O /tmp/intel-oneapi.gpg && \
    apt-key add /tmp/intel-oneapi.gpg && \
    rm -f /tmp/intel-oneapi.gpg && \
    echo "deb https://apt.repos.intel.com/oneapi all main" > /etc/apt/sources.list.d/oneAPI.list && \
    apt-get update && apt-get install -y --no-install-recommends \
    intel-oneapi-mkl \
    intel-oneapi-mkl-devel \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*
