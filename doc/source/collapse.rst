.. _collapse:


***************
Collapse fastq(.gz) files
***************

**Definition**

Normally quality values are lost in  small RNA-seq pipelines due to collapsing after adapter recognition. This option allow to collapse reads after adapter removal with ``cutadapt`` or any other tool. This way the mapping can use quality values, allowing to map using ``bwa`` for instance, or any other alignment tool that doesn't support FASTA files.

**Methods**

The new quality values are the average of each of the sequence collapse.

**Example**

::

    seqcluster collapse -f sample_trimmed.fastq -o collapse 

* ``-f`` is the fastq(.gz) file
* ``-o`` the folder where the outout will be created. A new FASTQ file, where the name stand for::

    @seq_[0-9]_x[0-9]


The number right after ``_x`` means the abundance of this sequence in the sample
