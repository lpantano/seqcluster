seqcluster
---------

Ask questions in the repo's Gitter: Join the chat at:

.. image:: https://badges.gitter.im/Join%20Chat.svg
    :target: https://gitter.im/lpantano/seqcluster
    
small RNA analysis from NGS data

.. image:: https://travis-ci.org/lpantano/seqcluster.png?branch=master
   :target: https://travis-ci.org/lpantano/seqcluster

.. image:: https://badge.fury.io/py/seqcluster.svg
   :target: http://badge.fury.io/py/seqcluster

.. image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square
   :target: http://bioconda.github.io
  
.. image:: https://anaconda.org/bioconda/seqcluster/badges/downloads.svg
   :target: https://anaconda.org/bioconda/seqcluster

.. image:: https://img.shields.io/pypi/l/seqcluster.svg
   :target: https://github.com/lpantano/seqcluster/blob/master/LICENSE
    
.. image:: http://www.repostatus.org/badges/latest/active.svg
   :alt: Project Status: Active - The project has reached a stable, usable state and is being actively developed.
   :target: http://www.repostatus.org/#active

.. image:: https://img.shields.io/badge/DOI-10.1093%2Fbioinformatics%2Fbtv632-lightgrey.svg?style=flat-square 
   :target: https://doi.org/10.1093/bioinformatics/btv632

Cite
---------

**Specific small-RNA signatures in the amygdala at premotor and motor stages of Parkinson's disease revealed by deep sequencing analysis**. Pantano L, Friedländer MR, Escaramís G, Lizano E, Pallarès-Albanell J, Ferrer I, Estivill X, Martí E.
Bioinformatics. 2015 Nov 2. pii: btv632. [Epub ahead of print]
PMID: `26530722 <http://www.ncbi.nlm.nih.gov/pubmed/26530722>`_

**A non-biased framework for the annotation and classification of the non-miRNA small RNA transcriptome.**
Pantano L1, Estivill X, Martí E. Bioinformatics. 2011 Nov 15;27(22):3202-3. doi: 10.1093/bioinformatics/btr527. Epub 2011 Oct 5.
PMID: `21976421 <http://www.ncbi.nlm.nih.gov/pubmed/21976421>`_

Quick start links
---------

See installation at http://seqcluster.readthedocs.org/installation.html

Moreover `bcbio-nextgen`_ provides 
a python framework to run a whole pipeline for small RNA (miRNA + tRNA + piRNA + others).

.. _bcbio-nextgen: https://bcbio-nextgen.readthedocs.org/en/latest/

An example of how to run with bcbio is here: http://seqcluster.readthedocs.org/example_pipeline.html#mirqc-data

In case you want to use seqcluster alone, a complete tutorial is here: http://seqcluster.readthedocs.org/getting_started.html#clustering-of-small-rna-sequences

Report
---------

Seqcluster creates html report that looks like `this`_. That is a table of all cluster detected, and you 
can go into each of them and get a complete `description`_ with profile expression figures, annotation details and
sequences counts for each sample. An interactive report is as well available to explore the expression profile,
position on the genome and secondary structure. Learn more http://seqcluster.readthedocs.org/more_outputs.html#report.

.. _this: https://rawgit.com/lpantano/seqcluster/master/data/examples_report/html/index.html
.. _description: https://rawgit.com/lpantano/seqcluster/master/data/examples_report/html/1/maps.html

Contributors
------------

* `Lorena Pantano  <https://github.com/lpantano>`_ (Bioinformatic Core, Harvard Chan School, Boston, USA)
* `Judith Flo Gaya <http://www.seas.harvard.edu/directory/jflo>`_ (School of Engineer and Aplied Science- Harvard University, Boston, USA)
* `Eulalia Marti Puig <http://www.crg.eu/en/group-members/eul%C3%A0lia-mart%C3%AD-puig>`_ (Genomics and Disease, Center of Genomic Regulation, Barcelona, Spain)
* Francisco Pantano Rubino: Architect

.. image:: https://d2weczhvl823v0.cloudfront.net/lpantano/seqcluster/trend.png
   :alt: Bitdeli badge
   :target: https://bitdeli.com/free

