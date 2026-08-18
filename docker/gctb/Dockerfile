# Builds the gctb binary FROM SOURCE against our GCTB fork, so we can patch
# upstream bugs (the make-ldm-eigen concurrency race, etc.). See
# plan_for_dummy_run_revised.txt and gctb_gwfm_constants.py for the fork/ref
# pin and the rationale. Used in genome-wide fine mapping with SBayesRC
# (gctb --gwfm RC); see https://gctbhub.cloud.edu.au/software/gctb/ .
#
# --------- build stage: compile gctb from the pinned fork commit ---------
FROM debian:stable-slim AS build
RUN apt-get update && apt-get install -y --no-install-recommends \
        g++ make git ca-certificates libeigen3-dev libboost-dev \
    && rm -rf /var/lib/apt/lists/*
ARG GCTB_FORK_REPO
ARG GCTB_FORK_REF
RUN git clone "$GCTB_FORK_REPO" /src \
    && git -C /src checkout "$GCTB_FORK_REF"
# Eigen3 and Boost are header-only (gctb links no compiled boost libs on Linux;
# LIB is just -ldl). The Makefile reads their include dirs from these vars.
RUN make -C /src/scr \
        EIGEN3_INCLUDE_DIR=/usr/include/eigen3 \
        BOOST_LIB=/usr/include \
    && test -x /src/scr/gctb
#
# --------- runtime stage: minimal image carrying just the binary ---------
FROM debian:stable-slim
# libstdc++6 + libgomp1 for the dynamically linked C++/OpenMP binary; unzip is
# invoked inside the container to extract the LD-matrix archives.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libstdc++6 libgomp1 ca-certificates unzip \
    && rm -rf /var/lib/apt/lists/*
COPY --from=build /src/scr/gctb /usr/local/bin/gctb
COPY --from=build /src/LICENSE /opt/gctb/LICENSE-GCTB
ENTRYPOINT []
