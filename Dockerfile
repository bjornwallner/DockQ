# DockQ container image (works locally and on RunPod pods)
#
# Build:  docker build -t dockq .
# Run:    docker run --rm -v "$PWD":/data dockq DockQ /data/model.pdb /data/native.pdb
#
# On RunPod the default command keeps the pod alive so you can open a
# terminal / SSH in and run `DockQ <model> <native>` directly.

FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

# Compiler is needed to build the Cython extension (DockQ.operations)
RUN apt-get update \
    && apt-get install -y --no-install-recommends build-essential \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/DockQ
COPY . .

RUN pip install --upgrade pip \
    && pip install . \
    && DockQ examples/1A2K_r_l_b.model.pdb examples/1A2K_r_l_b.pdb > /dev/null

# Working directory for user data (mount a RunPod volume here if desired)
WORKDIR /workspace

CMD ["sleep", "infinity"]
