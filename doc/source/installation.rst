.. _installation:

============
Installation
============

Seqcluster
----------

**With bcbio installed**

If you already have `bcbio`_, seqcluster comes with it. If you want the last development version::

/bcbio_anaconda_bin_path/seqcluster_install.py --upgrade

**Docker**::

    docker pull lpantano/smallsrna

**Bioconda binary**

install conda if you want an isolate env::

    wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh
    bash Miniconda-latest-Linux-x86_64.sh -b -p ~/install/seqcluster/anaconda


You can install directly from binstar (only for linux)::

    ~/install/seqcluster/anaconda/conda install seqcluster seqbuster bedtools samtools pip nose numpy scipy pandas pyvcf -c bioconda

With that you will have everything you need for the python package. 
The last step is to add seqcluster to your PATH if conda is not already there.

Go to Tools dependecies below to continue with the installation.

**Note**: After installation is highly recommended to get the last updated version doing::

    seqcluster_install.py --upgrade
   
**automated installation**

Strongly recommended to use `bcbio <https://bcbio-nextgen.readthedocs.org/en/latest/contents/installation.html>`_ installation if you work with sequencing data. But if you want a minimal installation::

    pip install fabric
    seqcluster_install --upgrade
    mkdir -p $PATH_TO_TOOLS/bin
    seqcluster_install --tools $PATH_TO_TOOLS

After that you will need to add to your path: ``export PATH=$PATH_TO_TOOLS/bin:$PATH``

Tools dependecies for a full small RNA pipeline
---------

For seqcluster command:

* bedtools
* samtools
* rnafold (for HTML report)

For some steps of a typical small RNA-seq pipeline (recommended to use directly `bcbio`_ ):

* STAR, bowtie
* fastqc
* cutadapt (install with ``bioconda`` using the same ``python`` env than seqcluster. 
You will need to link the ``cutadapt`` binary to your ``PATH``)
 
Data
---------

Easy way to install your small RNA seq data with `cloudbiolinux <https://github.com/chapmanb/cloudbiolinux>`_.
Seqcluster has snipped code to do that for you. Recommended to use `bcbio`_ for the pipeline since will install
everything you need in a single step ``bcbio_nextgen.py upgrade -u development --tools --genomes hg19 --aligners bowtie``.

But If you want to run ``seqcluster`` step by step an example of hg19 human version it will be (another well annotated supported genome is mm10):

Download genome data::

    seqcluster_install --data $PATH_TO_DATA --genomes hg19 --aligners bowtie2 --datatarget smallrna

If you want to install STAR indexes since gets kind of better results than bowtie2 (warning, 40GB memory RAM needed)::

    seqcluster_install --data $PATH_TO_DATA --genomes hg19 --aligners star


R package
---------

Install isomiRs package for R using devtools:: 

    devtools::install_github('lpantano/isomiRs')

To install all packages used by the Rmd report::

    Rscript -e 'source(https://raw.githubusercontent.com/lpantano/seqcluster/master/scripts/install_libraries.R)'


.. _bcbio: https://github.com/chapmanb/bcbio-nextgen
