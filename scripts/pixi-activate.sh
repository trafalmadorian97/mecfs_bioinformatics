# Sourced by pixi at activation time. Sets env vars that need an absolute path,
# which pixi's static [tool.pixi.activation].env table cannot express.

# R_LIBS_USER must be absolute so that R subprocesses spawned by
# parallel::makeCluster (PSOCK workers) re-evaluate it consistently. A relative
# path like "r-libs" would resolve relative to each worker's cwd, which can
# differ from the parent's, causing workers to fail to load packages installed
# under <project>/r-libs (e.g. GenomicSEM).
export R_LIBS_USER="${PIXI_PROJECT_ROOT}/r-libs"
