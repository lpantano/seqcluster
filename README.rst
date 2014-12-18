seqcluster
---------

small RNA analysis from NGS data

.. image:: https://travis-ci.org/lpantano/seqcluster.png?branch=master
    :target: https://travis-ci.org/lpantano/seqcluster.png?branch=master
    
.. image:: https://badge.fury.io/py/seqcluster.svg
    :target: http://badge.fury.io/py/seqcluster

.. image:: https://pypip.in/d/seqcluster/badge.png
    :target: https://pypi.python.org/pypi/seqcluster


installation
---------

Install first bcbio-nextgen and cutadapter::

    wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh
    bash Miniconda-latest-Linux-x86_64.sh -b -p ~/install/seqcluster/anaconda
    PATH = ~/install/seqcluster/anaconda/bin:$PATH
    conda install pip
    conda install -c https://conda.binstar.org/bcbio bcbio-nextgen
    pip install cutadapt

If you need to install bedtools, samtools and star::

   git clone https://github.com/Homebrew/linuxbrew.git  ~/install/seqcluster/linuxbrew
   cd ~/install/seqcluster/linuxbrew/bin
   ln -s `which gcc gcc-4.4`
   PATH = ~/install/seqcluster/linuxbrew/bin:$PATH
   brew tab homebrew/science
   brew tab chapmanb/homebre-cbl
   brew install bedtools
   brew install samtools
   brew install star
   

Then you can get seqcluster::

    pip install seqcluster

or the developement version::

    git clone https://github.com/lpantano/seqcluster
    cd seqcluster
    python setup.py install
