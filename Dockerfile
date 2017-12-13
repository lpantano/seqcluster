FROM ubuntu:latest
MAINTAINER Lorena Pantano "https://github.com/lpantano"


# Setup a base system 
RUN apt-get update && \
    apt-get install -y build-essential unzip wget git && \
    apt-get install -y libglu1-mesa && \
    apt-get install -y curl pigz bsdmainutils && \
    apt-get install -y --no-install-recommends libcurl4-gnutls-dev mbuffer python2.7-dev python-virtualenv 

RUN  wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh && \
     bash Miniconda-latest-Linux-x86_64.sh -b -p /usr/local/conda && \
     export PATH=/usr/local/conda/bin:$PATH && \
     conda config --add channels r && \
     conda install --yes ncurses -c r && \
     conda install --yes -c bioconda seqcluster seqbuster bedtools samtools pip nose numpy scipy pandas pyvcf cutadapt star -q

# setup paths
ENV  PATH="/usr/local/conda/bin:${PATH}"

