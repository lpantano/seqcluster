.. _installation:

============
Installation
============

Seqcluster
----------

**With bcbio installed**

If you already have `bcbio <https://github.com/chapmanb/bcbio-nextgen>`_, seqcluster comes with it. If you want the last development version just clone the repository.

**Binstar binary**

install conda if you want an isolate env::

    wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh
    bash Miniconda-latest-Linux-x86_64.sh -b -p ~/install/seqcluster/anaconda


You can install directly from binstar (only for linux)::

    ~/install/seqcluster/anaconda/conda install -c lpantano seqcluster -c bcbio

With that you will have everything you need for the python package. 
The last step is to add seqcluster to your PATH (see below).

Go to Tools dependecies below to continue with the installation.

**Step by step**

If you want to install step by step from a new conda environment::    

    ~/install/seqcluster/anaconda/bin/conda install pip
    ~/install/seqcluster/anaconda/bin/conda install -c https://conda.binstar.org/bcbio bcbio-nextgen

Remember to add the new python into your path every time you want to use seqcluster. 
If you already have `conda` in your system, just type::

    ~/install/seqcluster/anaconda/conda install -c https://conda.binstar.org/bcbio bcbio-nextgen

Then you can get seqcluster::

    ~/install/seqcluster/anaconda/conda install -c  https://conda.anaconda.org/lpantano seqcluster -c  https://conda.binstar.org/bcbio

Link binary to brew installation or to any folder is already in your path::

    ln -s ~/install/seqcluster/anaconda/bin/seqcluster* ~/install/seqcluster/linuxbrew/bin/.

**Note**: After installation is highly recommended to get the last updated version doing::

    seqcluster_install.py --upgrade

Tools dependecies
---------

For seqcluster command:

* bedtools
* samtools

For seqcluster-helper pipeline (this framework is deprecated, please use `bcbio`_ ):

* STAR
* fastqc
* cutadapt (install with ``pip`` using the same ``python`` env than seqcluster. 
You will need to link the ``cutadapt`` binary to your ``PATH``)::
    
    ~/install/seqcluster/anaconda/bin/pip install cutadapt

**easy installation**

Strongly recommend use `bcbio <https://bcbio-nextgen.readthedocs.org/en/latest/contents/installation.html>`_ installation if you work with sequencing data. But if you want a minimal installation::


    pip install fabric
    seqcluster_install --upgrade
    mkdir -p $PATH_TO_TOOLS/bin
    seqcluster_install --tools $PATH_TO_TOOLS


After that you will need to add to your path: ``export PATH=$PATH_TO_TOOLS/bin:$PATH``


**step by step**

To install dependencies using ``homebrew`` follow these steps::

   git clone https://github.com/Homebrew/linuxbrew.git  ~/install/seqcluster/linuxbrew
   cd ~/install/seqcluster/linuxbrew/bin
   ln -s `which gcc gcc-4.4`
   PATH = ~/install/seqcluster/linuxbrew/bin:$PATH
   brew tap homebrew/science
   brew tap chapmanb/homebrew-cbl
   brew install bedtools
   brew install samtools
   brew install star-rna
   brew install bowtie2

Data
---------

Easy way to install your small RNA seq data with `cloudbiolinux <https://github.com/chapmanb/cloudbiolinux>`_.
Seqcluster has snipped code to do that for you.

An exmaple of hg19 human version it will be (another well annotated supported genome is mm10):

Download genome data::

    seqcluster_install --data $PATH_TO_DATA --genomes hg19 --aligners bowtie2

If you want to install STAR indexes since gets kind of better results than bowtie2 (warning, 40GB memory RAM needed)::

    seqcluster_install --data $PATH_TO_DATA --genomes hg19 --aligners star


seqcluster-helper
---------

**Note: be aware that we moved to `bcbio`_ and seqcluster-helper is deprecated.**

`seqcluster-helper`_ provides 
a python framework to run a whole pipeline for small RNA (miRNA + others).

Assuming you installed seqcluster as mentioned before, clone this repository and type::

    python setup.py install
    ln -s ~/install/seqcluster/anaconda/bin/seqcluster-helper.py ~/install/seqcluster/linuxbrew/bin/.
    ln -s ~/install/seqcluster/anaconda/bin/seqcluster-installer.py ~/install/seqcluster/linuxbrew/bin/.

if you get problem with pythonpy: `pip install pythonpy`

**check installation**

::
    
    seqcluster-installer.py --check 

will tell you if all dependencies are installed and ready to use the framework


R pakcage
---------

Install isomiRs package for R using devtools:: 

    devtools::install_github('lpantano/isomiRs', ref='develop')

To install all packages used by the Rmd report::

    Rscript -e 'source(https://raw.githubusercontent.com/lpantano/seqcluster/master/scripts/install_libraries.R)'
    
    
.. _seqcluster-helper: https://github.com/lpantano/seqcluster-helper/blob/master/README.md
.. _bcbio: https://github.com/chapmanb/bcbio-nextgen
