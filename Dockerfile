FROM docker.io/library/ubuntu:24.04 AS builder
# The image needs no tests. To build them (-DNO_TESTS=OFF), also install
# libgtest-dev, or git and ca-certificates so that GoogleTest is downloaded.
ARG percolator_cmake_args="-DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DNO_TESTS=ON"

RUN apt-get update \
 && apt-get install -y --no-install-recommends build-essential cmake \
 && rm -rf /var/lib/apt/lists/*

COPY / /percolator

WORKDIR /build/percolator

RUN cmake ${percolator_cmake_args} /percolator \
 && make -j"$(nproc)" \
 && make install

FROM docker.io/library/ubuntu:24.04 AS runtime

# OpenMP runtime
RUN apt-get update \
 && apt-get install -y --no-install-recommends libgomp1 \
 && rm -rf /var/lib/apt/lists/*

COPY --from=builder /usr/bin/percolator /usr/bin/qvality /usr/bin/
