FROM docker.io/library/ubuntu:24.04 AS builder
ARG percolator_cmake_args="-DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DGOOGLE_TEST=0"

RUN apt-get update && apt-get install -y 

# Remove the system gtest libraries (if installed)
RUN apt-get remove --purge -y libgtest-dev && apt-get autoremove -y

RUN apt-get install -y build-essential g++ make cmake gawk git ca-certificates


COPY / /percolator

RUN mkdir -p /release /build
RUN mkdir -p /build/percolator;

WORKDIR /build/percolator

RUN cmake ${percolator_cmake_args} /percolator
RUN make
RUN make install

FROM docker.io/library/ubuntu:24.04 AS runtime

# OpenMP runtime
RUN apt-get update && apt-get install -y libgomp1

COPY --from=builder /usr/bin/percolator /usr/bin/percolator