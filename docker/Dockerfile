# Dockerfile for Normaliz

FROM ubuntu:xenial

RUN apt-get update -qq \
    && apt-get -qq install -y \
    build-essential m4 cmake \
    libboost-dev \
    libgmp-dev \
    git \
    wget \
    unzip \
    sudo

RUN adduser --quiet --shell /bin/bash --gecos "norm user,101,," --disabled-password norm \
    && adduser norm sudo \
    && chown -R norm:norm /home/norm \
    && echo 'norm ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

USER norm
ENV HOME /home/norm
ENV PATH ${HOME}/bin:${PATH}

WORKDIR /tmp

RUN cd /tmp

# ENV NUMBER_CORES=$(sudo cat /proc/cpuinfo | grep processor | wc -l)
ENV NUMBER_CORES = 4

RUN sudo wget http://cocoa.dima.unige.it/cocoalib/tgz/CoCoALib-0.99550.tgz \
    && sudo tar -xf CoCoALib-0.99550.tgz \
    && cd CoCoALib-0.99550 \
    && sudo ./configure --threadsafe-hack\
    && sudo make library -j${number_cores} \
    && cd /tmp

RUN git clone https://github.com/Normaliz/Normaliz.git && \
    cd Normaliz && \
    git checkout $1

RUN cd Normaliz && \
    mkdir BUILD && \
    cd BUILD && \
    COCOA_DIR=/tmp/CoCoALib-0.99550 cmake ../source && \
    make -j${number_cores} && \
    sudo make install
    sudo mkdir /home/norm/example
    sudo cp ../example/*.in /home/norm/example

RUN cd /tmp  \
    && sudo rm -rf Normaliz \
    && sudo rm -rf CoCoALib-0.99550 \
    && sudo rm -rf CoCoALib-config-root-2017-10-11 \
    && sudo rm CoCoALib-0.99550.tgz
