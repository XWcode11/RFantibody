FROM nvidia/cuda:11.8.0-cudnn8-devel-ubuntu22.04
ARG DEBIAN_FRONTEND=noninteractive

RUN sed -i 's/archive.ubuntu.com/mirrors.aliyun.com/g' /etc/apt/sources.list && \
    sed -i 's/security.ubuntu.com/mirrors.aliyun.com/g' /etc/apt/sources.list && \
    apt-get update && apt-get install -y software-properties-common && \
    DEBIAN_FRONTEND=noninteractive add-apt-repository -y ppa:deadsnakes/ppa && \
    apt-get install --no-install-recommends -y python3.10 python3-pip pipx vim make wget

RUN alias "python"="python3.10"

# Make a virtual env that we can safely install into

RUN python3 -m venv /opt/venv
# Enable venv
ENV PATH="/opt/venv/bin:$PATH"

RUN pip install poetry -i https://pypi.tuna.tsinghua.edu.cn/simple/

# Set the working directory to the user's home directory
WORKDIR /home

ENV http_proxy=http://127.0.0.1:7897
ENV https_proxy=http://127.0.0.1:7897

ENTRYPOINT ["/bin/bash"]