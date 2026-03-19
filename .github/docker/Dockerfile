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
    libopenblas-dev \
    libmetis-dev \
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

# Build and install sequential MUMPS with CMake package config
RUN git clone --depth=1 https://github.com/giavancini/mumps.git /tmp/mumps-src && \
    cmake -S /tmp/mumps-src -B /tmp/mumps-build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=/opt/mumps \
    -DBUILD_SHARED_LIBS=ON \
    -DBUILD_SINGLE=ON \
    -DBUILD_DOUBLE=ON \
    -DBUILD_COMPLEX=ON \
    -DBUILD_COMPLEX16=ON \
    -DMUMPS_parallel=OFF \
    -DMUMPS_openmp=ON \
    -DMUMPS_intsize64=ON \
    -DMUMPS_metis=ON \
    -DMUMPS_find_SCALAPACK=OFF \
    -DMUMPS_scalapack=OFF && \
    cmake --build /tmp/mumps-build -j"$(nproc)" && \
    cmake --install /tmp/mumps-build && \
    rm -rf /tmp/mumps-src /tmp/mumps-build && \
    echo "/opt/mumps/lib" > /etc/ld.so.conf.d/mumps.conf && \
    ldconfig
