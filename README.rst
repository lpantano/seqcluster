seqcluster
---------

Ask questions in the repo's Gitter: Join the chat at:

.. image:: https://badges.gitter.im/Join%20Chat.svg
    :target: https://gitter.im/lpantano/seqcluster
    
small RNA analysis from NGS data

.. image:: https://travis-ci.org/lpantano/seqcluster.png?branch=master
    :target: https://travis-ci.org/lpantano/seqcluster.png?branch=master

.. image:: https://badge.fury.io/py/seqcluster.svg
    :target: http://badge.fury.io/py/seqcluster

.. image:: https://pypip.in/d/seqcluster/badge.png
    :target: https://pypi.python.org/pypi/seqcluster


Cite
---------

A non-biased framework for the annotation and classification of the non-miRNA small RNA transcriptome.
Pantano L1, Estivill X, Martí E. Bioinformatics. 2011 Nov 15;27(22):3202-3. doi: 10.1093/bioinformatics/btr527. Epub 2011 Oct 5.
PMID: 21976421

installation
---------

`seqcluster-helper`_ provides 
a python framework to run a whole pipeline for small RNA (miRNA + others).

Install first bcbio-nextgen and cutadapter after install conda if you want a isolate env::

    wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh
    bash Miniconda-latest-Linux-x86_64.sh -b -p ~/install/seqcluster/anaconda
    PATH = ~/install/seqcluster/anaconda/bin:$PATH
    conda install pip
    conda install -c https://conda.binstar.org/bcbio bcbio-nextgen
    pip install cutadapt
    pip install matplotlib
    pip install -U cython


Remember to add the new python into your path every time you want to user seqcluster. 
If you already have `conda` in your system, just type::

    conda install -c https://conda.binstar.org/bcbio bcbio-nextgen

If you need to install bedtools, samtools and star, follow these steps::

   git clone https://github.com/Homebrew/linuxbrew.git  ~/install/seqcluster/linuxbrew
   cd ~/install/seqcluster/linuxbrew/bin
   ln -s `which gcc gcc-4.4`
   PATH = ~/install/seqcluster/linuxbrew/bin:$PATH
   brew tab homebrew/science
   brew tab chapmanb/homebre-cbl
   brew install bedtools
   brew install samtools
   brew install star-rna
   

Then you can get seqcluster::

    pip install seqcluster

or the developement version::

    git clone https://github.com/lpantano/seqcluster
    cd seqcluster
    python setup.py install


.. _seqcluster-helper: https://github.com/lpantano/seqcluster-helper/blob/master/README.md


quick start
---------

Complete tutorial here: http://seqcluster.readthedocs.org/getting_started.html#clustering-of-small-rna-sequences

report
---------

Seqcluster creates html report that looks like `this`_. That is a table of all cluster detected, and you 
can go into each of them and get a complete `description`_ with profile expression figures, annotation details and
sequences counts for each sample.

.. _this: https://rawgit.com/lpantano/seqcluster/master/data/examples_report/html/index.html
.. _description: https://rawgit.com/lpantano/seqcluster/master/data/examples_report/html/1/maps.html
