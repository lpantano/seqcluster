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
     conda install --yes -c conda-forge -c bioconda memory_profiler openjdk pysam pybedtools pandas numpy biopython progressbar2 pyyaml bedtools samtools mirtop viennarna pip nose setuptools -q && \
     wget https://github.com/lpantano/seqcluster/archive/master.zip && \
     unzip master.zip && \
     mv seqcluster-master seqcluster && \
     cd seqcluster && \
     /usr/local/conda/bin/python setup.py install && \
     /usr/local/conda/bin/nosetests
# setup paths
ENV  PATH="/usr/local/conda/bin:${PATH}"

