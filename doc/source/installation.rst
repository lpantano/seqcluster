.. _installation:


installation
--------

`seqcluster-helper`_ provides 
a python framework to run a whole pipeline for small RNA (miRNA + others).

Install first bcbio-nextgen and cutadapter after install conda if you want an isolate env::

    wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh
    bash Miniconda-latest-Linux-x86_64.sh -b -p ~/install/seqcluster/anaconda
    ~/install/seqcluster/anaconda/conda install pip
    ~/install/seqcluster/anaconda/conda install -c https://conda.binstar.org/bcbio bcbio-nextgen
    ~/install/seqcluster/anaconda/pip install cutadapt
    ~/install/seqcluster/anaconda/pip install matplotlib
    ~/install/seqcluster/anaconda/pip install -U cython


Remember to add the new python into your path every time you want to user seqcluster. 
If you already have `conda` in your system, just type::

    ~/install/seqcluster/anaconda/conda install -c https://conda.binstar.org/bcbio bcbio-nextgen

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

    ~/install/seqcluster/anaconda/pip install seqcluster

or the developement version::

    git clone https://github.com/lpantano/seqcluster
    cd seqcluster
    ~/install/seqcluster/anaconda/python setup.py install

Link binary to brew installation or to any folder is already in your path::

    ln -s ~/install/seqcluster/anaconda/bin/seqcluster ~/install/seqcluster/linuxbrew/bin/.

You can install the python framework for the full small RNA analysis (`seqcluster-helper`_)::

    brew install https://github.com/lpantano/seqcluster-helper/blob/master/seqbuster.rb
    brew install fastqc

And finally clone this repository and type::

    python setup.py install
    ln -s ~/install/seqcluster/anaconda/bin/seqcluster-helper.py ~/install/seqcluster/linuxbrew/bin/.
    ln -s ~/install/seqcluster/anaconda/bin/seqcluster-installer.py ~/install/seqcluster/linuxbrew/bin/.


if you get problem with pythonpy: `pip install pythonpy`

Install isomiRs package for R using devtools:: 

    devtools::install_github('lpantano/isomiRs', ref='develop')


.. _seqcluster-helper: https://github.com/lpantano/seqcluster-helper/blob/master/README.md


check installation
--------

`seqcluster-installer.py --check` will tell you if all dependencies are installed and ready to use the framework

Easy star seqcluster-helper.py
--------


`seqcluster-helper.py --sample-map config.csv --aligner-index /path/2/star_index --gtf-file /path/2/gtf_annotation --species hsa --reference /path/2/genome/genome.fasta`

* `sample-map` file should be a csv file with: `name,/path/2/fastq,group` for each sample
* `genome.fasta` needs to have the FAI file. You can create this with: `samtools faidx genome.fasta`
* `gtf-file` is used for annotation. The 3 column is the group of sRNA and the `gene_name` attribute the annotation
* `species` should be compatible with miRBase notation
* `DB` is the path to `harpin.fa` and `miRNAstr`, like this https://github.com/lpantano/seqbuster/tree/master/modules/miraligner/DB

Options to run in a cluster
--------

It uses ipython-cluster-helper to send jobs to nodes in the cluster

* `--parallel` should set to `ipython`
* `--scheduler` should be set to `sge,lsf,slurm`
* `--num-jobs` indicates how much jobs to launch. It will run samples independently. If you have 4 samples, and set this to 4, 4 jobs will be launch to the cluster
* `--queue` the queue to use
* `--resources` allows to set any special parameter for the cluster, such as, email in sge system: `M=my@email.com`

Read complete usability here: https://github.com/roryk/ipython-cluster-helper
An examples in slurm system is `--parallel ipython --scheduler slurm --num-jobs 4 --queue general`

Output
--------

* one folder for each sample
 * adapter: `*clean.fastq` is the file after adapter removal, `*clean_trimmed.fastq` is the collapse `clean.fastq`, `*fragments.fastq` is file without adapter, `*short.fastq` is file with reads < 16 nt.
 * align: BAM file results from align `trimmed.fastq`
 * miraligner: file with miRNA anotation 
 * qc: `*_fastqc.html` is the fastqc results from the uncollapse fastq file
* seqcluster: is the result of running seqcluster. See its [documentation](http://seqcluster.readthedocs.org/getting_started.html#clustering-of-small-rna-sequences) for further information.
* `report-ready.Rmd`: template to create a quick html report with exploration and differential expression analysis.
