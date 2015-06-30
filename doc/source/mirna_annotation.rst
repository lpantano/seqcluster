.. _mirna_annotation:


***************
miRNA annotation
***************


**REMOVE ADAPTER**

I am currently using ``cutadapt``.

::

    cutadapt --adapter=$ADAPTER --minimum-length=8 --untrimmed-output=sample1_notfound.fastq -o sample1_clean.fastq -m 17 --overlap=8 sample1.fastq 

**COLLAPSE READS**

To reduce computational time, I recommend to collapse sequences, also it would help to apply filters based on abundances.
Like removing sequences that appear only once.

::

   seqcluster collapse -f sample1_clean.fastq -o collapse

Here I am only using sequences that had the adapter, meaning that for sure are small fragments.

**MIRALIGNER**

Download the tool from `miraligner`_ repository. 

.. _miraligner: https://github.com/lpantano/seqbuster/blob/master/modules/miraligner/miraligner.jar

Download the mirbase files (`hairpin`_ and `miRNA`_) from the ftp and save it to `DB` folder.

.. _hairpin: ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.zip
.. _miRNA: ftp://mirbase.org/pub/mirbase/CURRENT/miRNA.str.zip

You can map the miRNAs with.

::

     java -jar miraligner.jar -sub 1 -trim 3 -add 3 -s hsa -i test.fa -db DB  -o output_prefix 

You can read more about the output here in the `wiki`_

.. _wiki: https://github.com/lpantano/seqbuster/wiki/miraligner
