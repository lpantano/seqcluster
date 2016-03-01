FROM stackbrew/ubuntu:14.04
MAINTAINER Lorena Pantano "https://github.com/lpantano"

# Release 0.9.5a -- https://github.com/chapmanb/bcbio-nextgen/commit/40fb5ba

# Setup a base system 
RUN apt-get update && apt-get install -y build-essential zlib1g-dev wget curl python-setuptools git && \
    apt-get install -y openjdk-7-jdk openjdk-7-jre ruby libncurses5-dev libcurl4-openssl-dev libbz2-dev \
    unzip pigz bsdmainutils leiningen && \

# Fake a fuse install; openjdk pulls this in 
# https://github.com/dotcloud/docker/issues/514
# https://gist.github.com/henrik-muehe/6155333
    mkdir -p /tmp/fuse-hack && cd /tmp/fuse-hack && \
    apt-get install libfuse2 && \
    apt-get download fuse && \
    dpkg-deb -x fuse_* . && \
    dpkg-deb -e fuse_* && \
    rm fuse_*.deb && \
    echo -en '#!/bin/bash\nexit 0\n' > DEBIAN/postinst && \
    dpkg-deb -b . /fuse.deb && \
    dpkg -i /fuse.deb && \
    rm -rf /tmp/fuse-hack && \

# bcbio-nextgen installation
    mkdir -p /tmp/bcbio-nextgen-install && cd /tmp/bcbio-nextgen-install && \
    wget --no-check-certificate \
      https://raw.github.com/chapmanb/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py && \
    python bcbio_nextgen_install.py /usr/local/bcbio \
      --nodata -u development --isolate && \

# setup paths
    echo 'export PATH=/usr/local/bcbio/anaconda/bin:$PATH' >> /etc/profile.d/bcbio.sh && \

    /usr/local/bcbio/anaconda/bin/seqcluster_install --upgrade
    mkdir -p /usr/local/tools/bin
    /usr/local/bcbio/anaconda/bin/seqcluster_install --tools /usr/local/tools
# add user run script
    wget --no-check-certificate -O createsetuser \
      https://raw.github.com/chapmanb/bcbio-nextgen-vm/master/scripts/createsetuser && \
    chmod a+x createsetuser && mv createsetuser /sbin && \

# clean filesystem
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/tmp/* && \
    /usr/local/bcbio/anaconda/bin/conda remove --yes qt && \
    /usr/local/bcbio/anaconda/bin/conda clean --yes --tarballs && \
    rm -rf /usr/local/anaconda/pkgs/qt* && \
    rm -rf $(brew --cache) && \
    rm -rf /.cpanm && \
    rm -rf /tmp/bcbio-nextgen-install && \

RUN export PATH=/usr/local/bcbio/anaconda/bin:/usr/local/tools:$PATH
# if you want to install data:
# seqcluster_install --data /usr/local/bcbio --genomes hg19 --aligners bowtie2
