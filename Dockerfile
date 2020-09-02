FROM ubuntu:latest

RUN apt-get -y update 

RUN apt-get install -y sudo wget make git curl liblzma-dev libz-dev libncurses5-dev libncursesw5-dev libghc-bzlib-prof gcc unzip zlib1g-dev libncurses5-dev libbz2-dev liblzma-dev libcurl4-openssl-dev python3 python3-numpy python3-pip

RUN pip3 install biopython

COPY ./Install_dependancies.sh /Install_dependancies.sh

# RUN useradd -m docker && echo "docker:docker" | chpasswd && adduser docker sudo
# USER docker

RUN adduser --disabled-password --gecos '' docker
RUN adduser docker sudo
RUN echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

USER docker

RUN sudo ln -s /usr/bin/python3 /usr/bin/python
RUN sudo mkdir -p /BovTB-nf/

RUN sudo sh /Install_dependancies.sh