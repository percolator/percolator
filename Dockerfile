FROM docker.io/library/ubuntu:24.04 as builder
ARG percolator_cmake_args="-DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DGOOGLE_TEST=0 -DXML_SUPPORT=OFF"

RUN apt-get update && apt-get install -y 

RUN apt-get install -y g++ make cmake gawk

RUN apt-get install -y \
      libboost-system1.74.0 \
      libboost-thread1.74.0 

RUN apt-get install -y -o Acquire::Retries=3 libxml2-utils

RUN mkdir -p /release /build
COPY / /percolator

RUN mkdir -p /build/percolator-noxml /build/percolator /build/converters;

WORKDIR /build/percolator-noxml

RUN cmake ${percolator_cmake_args} /percolator
RUN make
RUN make install

FROM docker.io/library/ubuntu:24.04

RUN apt-get update && apt-get install -y \
      libboost-system1.74.0 \
      libboost-thread1.74.0 \
      libbz2-1.0 \
      libcurl4t64 \
      libgomp1 \
      zlib1g

COPY --from=builder /usr/bin/percolator /usr/bin/percolator