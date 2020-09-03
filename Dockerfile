FROM ubuntu:latest

# apt-get dependencies
RUN apt-get -y update 
RUN apt-get install -y sudo \
    wget \
    make \
    git \
    curl \
    liblzma-dev \
    libz-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libghc-bzlib-prof \
    gcc unzip zlib1g-dev \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    python3 \
    python3-numpy \
    python3-pip

# Setup Sudo. Thanks: https://stackoverflow.com/questions/25845538/how-to-use-sudo-inside-a-docker-container
RUN adduser --disabled-password --gecos '' docker
RUN adduser docker sudo
RUN echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers
USER docker

# Setup python 
RUN pip3 install biopython
RUN sudo ln -s /usr/bin/python3 /usr/bin/python

# Install Dependencies
RUN sudo mkdir -p /BovTB-nf/
COPY ./Install_dependancies.sh /Install_dependancies.sh
RUN sudo sh /Install_dependancies.sh
