FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install necessary packages for building and running the application
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    g++ \
    gfortran \
    cmake \
    ninja-build \
    git \
    liblapack-dev \
    libblas-dev \
    libboost-all-dev \
    python3 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*