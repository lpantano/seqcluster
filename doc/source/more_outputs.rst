.. _more_output:


***************
More outputs
***************

You can obtain different statistics from the analyais:

* abundance distribution by length of reads that have been aligned, 
that have appear in the output, and annotated with different databases.
To obtain this, you need to run this command::

    samtools index Aligned.sortedByCoord.out.bam
    seqcluster stats -j res/cluster/seqcluster.json -m res/seqs.ma -a Aligned.sortedByCoord.out.bam -o res/stats 
    
The output file is ``stats_align.dat``. It is a 4 column file with the following information:

* size of the read
* sample
* expression
* class: ALIGN (aligned  read), JSOB (included in final output), DATABASE (names of the databases assigned to the read)
