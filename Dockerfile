FROM docker.io/library/ubuntu:24.04 as builder
ARG percolator_cmake_args="-DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DXML_SUPPORT=OFF"

RUN apt-get update && apt-get install -y 

RUN apt-get install -y -o Acquire::Retries=3 libxml2-utils
RUN mkdir -p /release /build
COPY / /percolator

WORKDIR /

RUN ./admin/builders/ubuntu64_build.sh -s / -r /release -b /build

FROM docker.io/library/ubuntu:24.04

RUN apt-get update && apt-get install -y \
      libboost-system1.74.0 \
      libboost-thread1.74.0 \
      libbz2-1.0 \
      libcurl4t64 \
      libgomp1 \
      libsqlite3-0 \
      zlib1g
COPY --from=builder /build/percolator-noxml/src/percolator /usr/bin/percolator