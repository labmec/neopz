# =============================================================
# NeoPZ CI dependencies image  (neopz-deps)
#
# Provides all build dependencies pre-installed so that CI jobs
# only need to checkout + cmake + ctest without spending time
# downloading/building libraries on every run.
#
# NeoPZ source is NOT included — GitHub Actions mounts the
# checked-out workspace into the container automatically.
#
# Always installed:
#   - Build tools: cmake, ninja, gcc, g++, gfortran, git, etc.
#   - libmetis-dev, libboost-all-dev, liblog4cxx-dev
#   - Intel oneAPI MKL / Pardiso  (x86_64)
#   - MUMPS (d+z always; s and c via build args)
#
# Build args (0 = off, 1 = on):
#   MUMPS_BUILD_SINGLE    – include MUMPS s (real float) variant  (default: 0)
#   MUMPS_BUILD_COMPLEX   – include MUMPS c (complex float) variant (default: 0)
#   MUMPS_BUILD_COMPLEX16 – include MUMPS z (complex double) variant (default: 1)
#
# Environment variables set in the final image:
#   MUMPS_ROOT – /opt/mumps
#   MKLROOT    – /opt/intel/oneapi/mkl/latest
# =============================================================


# =============================================================
# Stage 1 – MUMPS builder
# Builds MUMPS from the giavancini/mumps CMake-based fork.
# Variant flag names follow MUMPS/LAPACK conventions:
#   BUILD_DOUBLE=on    → dmumps (d - real double)      always ON
#   BUILD_SINGLE=on    → smumps (s - real single)
#   BUILD_COMPLEX=on   → cmumps (c - complex float)
#   BUILD_COMPLEX16=on → zmumps (z - complex double)
# =============================================================
FROM debian:trixie-slim AS builder-mumps

ARG MUMPS_BUILD_SINGLE=0
ARG MUMPS_BUILD_COMPLEX=0
ARG MUMPS_BUILD_COMPLEX16=1
ARG MUMPS_REPO=https://github.com/giavancini/mumps.git

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential cmake git gfortran \
        libopenblas-dev libmetis-dev \
    && rm -rf /var/lib/apt/lists/*

RUN git clone --depth=1 "${MUMPS_REPO}" /tmp/mumps-src && \
    cmake -S /tmp/mumps-src -B /tmp/mumps-src/build \
        -DBUILD_DOUBLE=on \
        "-DBUILD_SINGLE=$([ "${MUMPS_BUILD_SINGLE}"      = "1" ] && echo on || echo off)" \
        "-DBUILD_COMPLEX=$([ "${MUMPS_BUILD_COMPLEX}"     = "1" ] && echo on || echo off)" \
        "-DBUILD_COMPLEX16=$([ "${MUMPS_BUILD_COMPLEX16}" = "1" ] && echo on || echo off)" \
        -DMUMPS_parallel=false \
        -DMUMPS_openmp=on \
        -DBUILD_SHARED_LIBS=on \
        -DMUMPS_intsize64=on \
        -DMUMPS_metis=on \
        -DMUMPS_find_SCALAPACK=false \
        -DMUMPS_scalapack=false \
        -DLAPACK_VENDOR=OpenBLAS \
        -DCMAKE_INSTALL_PREFIX=/opt/mumps && \
    cmake --build /tmp/mumps-src/build -j$(nproc) && \
    cmake --install /tmp/mumps-src/build && \
    rm -rf /tmp/mumps-src


# =============================================================
# Stage 2 – Final CI image
# All dependencies ready; NeoPZ is compiled by the CI job.
# =============================================================
FROM debian:trixie-slim

ENV DEBIAN_FRONTEND=noninteractive

# ── Core build tools + always-on optional libraries ──────────
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        g++ \
        gfortran \
        cmake \
        ninja-build \
        git \
        wget \
        ca-certificates \
        gnupg \
        pkg-config \
        libopenblas-dev \
        liblapack-dev \
        libmetis-dev \
        libboost-all-dev \
        liblog4cxx-dev \
        python3 \
    && rm -rf /var/lib/apt/lists/*

# ── Intel oneAPI MKL / Pardiso ────────────────────────────────
RUN wget -qO- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
        | gpg --dearmor -o /usr/share/keyrings/intel-oneapi-archive-keyring.gpg && \
    echo "deb [signed-by=/usr/share/keyrings/intel-oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" \
        > /etc/apt/sources.list.d/oneAPI.list && \
    apt-get update && \
    apt-get install -y --no-install-recommends intel-oneapi-mkl intel-oneapi-mkl-devel && \
    rm -rf /var/lib/apt/lists/*

# ── Import MUMPS artifacts from builder stage ─────────────────
COPY --from=builder-mumps /opt/mumps /opt/mumps

# ── Register shared libraries ─────────────────────────────────
RUN echo "/opt/mumps/lib" > /etc/ld.so.conf.d/mumps.conf && \
    echo "/opt/intel/oneapi/mkl/latest/lib/intel64" > /etc/ld.so.conf.d/mkl.conf && \
    ldconfig

# ── Persistent environment variables for cmake find_package ──
ENV MUMPS_ROOT=/opt/mumps
ENV MKLROOT=/opt/intel/oneapi/mkl/latest

WORKDIR /workspace
CMD ["/bin/bash"]