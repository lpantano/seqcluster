.. _installation:

============
Installation
============

Seqcluster
----------

**With bcbio installed**

If you already have `bcbio`_, seqcluster comes with it. If you want the last development version::

/bcbio_anaconda_bin_path/seqcluster_install.py --upgrade

**Binstar binary**

install conda if you want an isolate env::

    wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh
    bash Miniconda-latest-Linux-x86_64.sh -b -p ~/install/seqcluster/anaconda


You can install directly from binstar (only for linux)::

    ~/install/seqcluster/anaconda/conda install seqcluster bcbio-nextgen -c bioconda -c bcbio

With that you will have everything you need for the python package. 
The last step is to add seqcluster to your PATH (see below).

Go to Tools dependecies below to continue with the installation.

**Step by step**

If you want to install step by step from a new conda environment::    

    ~/install/seqcluster/anaconda/bin/conda install pip
    ~/install/seqcluster/anaconda/bin/conda install -c bcbio bcbio-nextgen

Link binary to any folder is already in your path::

    ln -s ~/install/seqcluster/anaconda/bin/seqcluster* ~/install/seqcluster/linuxbrew/bin/.

**Note**: After installation is highly recommended to get the last updated version doing::

    seqcluster_install.py --upgrade

Tools dependecies
---------

For seqcluster command:

* bedtools
* samtools

For some steps of a typical small RNA-seq pipeline (recommended to use directly `bcbio`_ ):

* STAR, bowtie
* fastqc
* cutadapt (install with ``pip`` using the same ``python`` env than seqcluster. 
You will need to link the ``cutadapt`` binary to your ``PATH``)
    
**easy installation**

Strongly recommended to use `bcbio <https://bcbio-nextgen.readthedocs.org/en/latest/contents/installation.html>`_ installation if you work with sequencing data. But if you want a minimal installation::

    pip install fabric
    seqcluster_install --upgrade
    mkdir -p $PATH_TO_TOOLS/bin
    seqcluster_install --tools $PATH_TO_TOOLS

After that you will need to add to your path: ``export PATH=$PATH_TO_TOOLS/bin:$PATH``


Data
---------

Easy way to install your small RNA seq data with `cloudbiolinux <https://github.com/chapmanb/cloudbiolinux>`_.
Seqcluster has snipped code to do that for you. Recommended to use `bcbio`_ for the pipeline since will install
everything you need in a single step ``bcbio_nextgen.py upgrade -u development --tools --genomes hg19 --aligners bowtie``.

But If you want to run ``seqcluster`` step by step an exmaple of hg19 human version it will be (another well annotated supported genome is mm10):

Download genome data::

    seqcluster_install --data $PATH_TO_DATA --genomes hg19 --aligners bowtie2

If you want to install STAR indexes since gets kind of better results than bowtie2 (warning, 40GB memory RAM needed)::

    seqcluster_install --data $PATH_TO_DATA --genomes hg19 --aligners star


R pakcage
---------

Install isomiRs package for R using devtools:: 

    devtools::install_github('lpantano/isomiRs')

To install all packages used by the Rmd report::

    Rscript -e 'source(https://raw.githubusercontent.com/lpantano/seqcluster/master/scripts/install_libraries.R)'


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

    
    
.. _seqcluster-helper: https://github.com/lpantano/seqcluster-helper/blob/master/README.md
.. _bcbio: https://github.com/chapmanb/bcbio-nextgen
