# This dockerfile is used in genome-wide fine mapping with SBayesRC using the gctb software
# See https://gctbhub.cloud.edu.au/software/gctb/ for information on the software
FROM debian:stable-slim
RUN apt-get update && apt-get install -y --no-install-recommends \
        libgomp1 ca-certificates curl unzip \
    && rm -rf /var/lib/apt/lists/*
ARG GCTB_URL
ARG GCTB_SHA256
RUN curl -fsSL "$GCTB_URL" -o /tmp/gctb.zip \
    && echo "$GCTB_SHA256  /tmp/gctb.zip" | sha256sum -c - \
    && unzip /tmp/gctb.zip -d /opt/gctb && rm /tmp/gctb.zip \
    && ln -s "$(find /opt/gctb -name gctb -type f | head -1)" /usr/local/bin/gctb
COPY LICENSE-GCTB /opt/gctb/LICENSE-GCTB
ENTRYPOINT []
