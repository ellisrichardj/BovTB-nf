FROM ubuntu:20.04
WORKDIR /BovTB-nf/

################## METADATA ##########################

LABEL base.image="ubuntu:20.04"
LABEL software="Bovine-TB Pipeline Image"
LABEL about.summary="Bioinformatics Pipeline for post-processing of Bovine TB fastq reads"
LABEL about.documentation="https://github.com/ellisrichardj/BovTB-nf"
LABEL about.tags="Genomics, WGS"


################## DEPENDENCIES ######################

# apt-get dependencies
RUN apt-get -y update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    openjdk-8-jdk \
    sudo \
    wget \
    make \
    git \
    curl \
    liblzma-dev \
    libz-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libghc-bzlib-prof \
    gcc \
    unzip \
    zlib1g-dev \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    python3 \
    python3-numpy \
    python3-pip

# Install nextflow.
# ENV line required to set jvm memory and cpu limits to docker. 
# see: https://github.com/nextflow-io/nextflow/blob/v20.07.1/docker/Dockerfile
ENV NXF_OPTS='-XX:+UnlockExperimentalVMOptions -XX:+UseCGroupMemoryLimitForHeap' NXF_HOME=/.nextflow
COPY install_nextflow-20.7.1.bash ./install_nextflow.bash
RUN cat ./install_nextflow.bash | bash
RUN ln -s $PWD/nextflow /usr/local/bin/nextflow

# python 
RUN pip3 install biopython
RUN ln -s /usr/bin/python3 /usr/bin/python

# bovtb tools
COPY ./Install_dependancies.sh ./Install_dependancies.sh
RUN sh ./Install_dependancies.sh

# pipeline
COPY ./bTB-WGS_process.nf ./


################## ENTRY ######################

CMD /bin/bash
