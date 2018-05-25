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
     conda install --yes -c conda-forge -c bioconda seqcluster bedtools samtools pip nose numpy scipy pandas pyvcf pytz dateutil setuptools -q
RUN wget https://github.com/lpantano/seqcluster/archive/master.zip && \
    unzip master.zip && \
    cd seqcluster-master && \
    /usr/local/conda/bin/python setup.py install
# setup paths
ENV  PATH="/usr/local/conda/bin:${PATH}"

