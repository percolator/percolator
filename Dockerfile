FROM docker.io/library/ubuntu:24.04 as builder
ARG percolator_cmake_args="-DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DXML_SUPPORT=OFF"

RUN apt-get update && apt-get install -y 

RUN apt-get install -y -o Acquire::Retries=3 libxml2-utils

# Create a "fake sudo" script that simply strips off the word 'sudo'
# and then runs the rest of the command as-is.
RUN echo '#!/bin/bash\nexec "$@"' > /usr/bin/sudo && \
    chmod +x /usr/bin/sudo

RUN apt-get install -y g++ make cmake gawk

RUN mkdir -p /release /build
COPY / /percolator

WORKDIR /

RUN /percolator/admin/builders/ubuntu64_build.sh -s / -r /release -b /build

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