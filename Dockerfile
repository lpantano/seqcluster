FROM bcbio/bcbio:latest
MAINTAINER Lorena Pantano "https://github.com/lpantano"

# Setup a base system 
RUN /usr/local/share/bcbio-nextgen/anaconda/bin/bcbio_nextgen.py upgrade -u development --tools && \ 
    /usr/local/share/bcbio-nextgen/anaconda/bin/seqcluster_install --upgrade &&\
    ln -s /usr/local/share/bcbio-nextgen/anaconda/bin/seqcluster* /usr/local/bin/.
# if you want to install data:
# seqcluster_install --data /usr/local/bcbio --genomes hg19 --aligners bowtie2
