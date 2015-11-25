.. _mirna_annotation:


***************
miRNA annotation
***************

miRNA annotation is running inside `bcbio small RNAseq pipeline <https://bcbio-nextgen.readthedocs.org/en/latest/contents/pipelines.html#smallrna-seq>`_ together with other tools to do a complete
small RNA analysis.

For some comparison with other tools go `here <https://github.com/lpantano/mypubs/blob/master/mirna/mirannotation/stats.md>`_.

You can run samples after processing the reads as shown below.
Currently there are two version: JAVA and PYTHON/C. 

**Naming**

It is a working process, but since 10-21-2015 isomiR naming has changed to:

* Changes at 5' end: ``0/NA`` means no modification. ``UPPER CASE LETTER`` means nucleotide insertions (sequence starts before miRBase mature position). ``LOWWER CASE LETTER`` means nucleotide deletions (sequence starts after miRBase mature position).
* Changes at 3' end: ``0/NA`` means no modification. ``UPPER CASE LETTER`` means nucleotide insertions (sequence ends after miRBase mature position). ``LOWWER CASE LETTER`` means nucleotide deletions (sequence ends before miRBase mature position).
* Additions at 3' end: ``0/NA`` means no modification. ``UPPER CASE LETTER`` means addition at the end. Note these nucleotides don't match the precursor. So they are post-transcriptional modification.
* Nucleotide substitution: ``NUMBER|NUCLEOTIDE_ISOMIR|NUCLEOTIDE_REFERENCE`` means at the position giving by the number the nucleotide in the sequence has substituted the nucleotide in the reference. This, as well, is a post-transcriptional modification.

There are two different outputs right now: 1) tab limited format, where each column will refer to the previous 4 points, or 2) a merged format, where these 4 points (mirna, substitution, addition, 5 trimming, 3 trimming) are separated by ``:``. For instance: ``hsa-let-7a-5p:0:0:GT:t`` means ``hsa-let-7a-5p`` has a 5' trimming event (tarts 2 nts before the reference miRNA) and a 3' trimming event (ends 1 nt before the reference miRNA).

Processing of reads
-------------------

**REMOVE ADAPTER**

I am currently using ``cutadapt``.

::

    cutadapt --adapter=$ADAPTER --minimum-length=8 --untrimmed-output=sample1_notfound.fastq -o sample1_clean.fastq -m 17 --overlap=8 sample1.fastq 

**COLLAPSE READS**

To reduce computational time, I recommend to collapse sequences, also it would help to apply filters based on abundances.
Like removing sequences that appear only once.

::

   seqcluster collapse -f sample1_clean.fastq -o collapse

Here I am only using sequences that had the adapter, meaning that for sure are small fragments. The output will be named as ``sample1_clean_trimmed.fastq``


Prepare databases
-----------------

For human or mouse, follows `this instruction <http://seqcluster.readthedocs.org/installation.html#data>`_ to download easily miRBase files. For other species you only need hairpin.fa and miRNA.str from miRBase site. **Highly recommended to filer hairpin.hsa to contain only the desired species.**

miRNA/isomiR annotation with JAVA
---------------------------------

**MIRALIGNER**

Download the tool from `miraligner`_ repository. 

.. _miraligner: https://github.com/lpantano/seqbuster/blob/master/modules/miraligner/miraligner.jar

Download the mirbase files (`hairpin`_ and `miRNA`_) from the ftp and save it to `DB` folder.

.. _hairpin: ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.zip
.. _miRNA: ftp://mirbase.org/pub/mirbase/CURRENT/miRNA.str.zip

You can map the miRNAs with.

::

     java -jar miraligner.jar -sub 1 -trim 3 -add 3 -s hsa -i sample1_clean_trimmed.fastq -db DB  -o output_prefix 


**Cite**

SeqBuster is a bioinformatic tool for the processing and analysis of small RNAs datasets, reveals ubiquitous miRNA modifications in human embryonic cells. Pantano L, Estivill X, Martí E. *Nucleic Acids Res. 2010 Mar;38(5):e34. Epub 2009 Dec 11.*

miRNA/isomiRs annotation with python
------------------------------------

A new function to annotate miRNA/isomiR sequences using BAM files aligned to miRBase precursors or fastq files to align from scratch has been added to ``seqcluster``::

	seqcluster seqbuster --out results --hairpin hairpin.fa --mirna miRNA.str --species hsa input_file.fastq ...

If the input file is a BAM file, seqcluster will parse it to produce miRNA annotation, including isomiRs. If the input is FASTQ/FASTA
file, seqcluster will map with the new C implementation of miraligner and annotate miRNAs and isomiRs as before. 

Multiple files can be given to analyze all of them serially. Files inside the output folder are:

* raw mirna annotation to all posible mirnas (``*.premirna``) 
* count file for miRNAs (``counts_mirna.tsv``) 
* count file for isomiRs (``counts.tsv``) 

**NOTE:** `Check comparison of multiple tools <https://github.com/lpantano/mypubs/blob/master/mirna/mirannotation/stats.md>`_ for miRNA annotation.

Post-analysis with R
--------------------

Use the outputs to do differential expression, clustering and descriptive analysis with this package: `isomiRs <https://github.com/lpantano/isomiRs>`_

Manual of miraligner(JAVA)
--------------------------

**options**

Add ``-freq`` if you have your fasta/fastq file with this format and you want a third column with the frequency (normally value after x character)::


    >seq_1_x4
    CACCGCTGTCGGGGAACCGCGCCAATTT


Add ``-pre`` if you want also sequences that map to the precursor but outside the mature miRNA


* Parameter `-sub`: mismatches allowed (0/1)
* Parameter `-trim`: nucleotides allowed for trimming (max 3)
* Parameter `-add`: nucleotides allowed for addition (max 3)
* Parameter `-s`: species (3 letter, human=>hsa)
* Parameter `-i`: fasta file
* Parameter `-db`: folder where miRBase files are(one copy at miraligner-1.0/DB folder)
* Parameter `-o`: prefix for the output files
* Parameter `-freq`: add frequency of the sequence to the output (just where input is fasta file with name matching this patter: >seq_3_x67)
* Parameter `-pre`: add sequences mapping to precursors as well

**input**

A fasta/fastq file reads::

    >seq
    CACCGCTGTCGGGGAACCGCGCCAATTT

or tabular file with counts information::

CACCGCTGTCGGGGAACCGCGCCAATTT 45

**output**

Track file *.mirna.opt: information about the process

Non mapped sequences will be on *.nomap

Header of the *.mirna.out file:

* seq: sequence
* freq/name: depending on the input this column contains counts (tabular input file) or name (fasta file)
* mir: miRNA name
* start: start of the sequence at the precursor
* end: end of the sequence at the precursor
* mism: nucleotide substitution position | nucleotide at sequence | nucleotide at precursor
* addition: nucleotides at 3 end added::


    precursor         => cctgtggttagctggttgcatatcc
    annotated miRNA   =>   TGTGGTTAGCTGGTTGCATAT
    sequence add:  TT =>   TGTGGTTAGCTGGTTGCATATTT


* tr5: nucleotides at 5 end different from the annonated sequence in miRBase::


	precursor 	      => cctgtggttagctggttgcatatcc
	annotated miRNA   =>   TGTGGTTAGCTGGTTGCATAT
	sequence tr5:  CC => CCTGTGGTTAGCTGGTTGCATAT
	sequence tr5:  tg =>     TGGTTAGCTGGTTGCATAT


* tr3: nucleotides at 3 end different from the annotated sequence in miRBase::


    precursor         => cctgtggttagctggttgcatatcc
    annotated miRNA   =>   TGTGGTTAGCTGGTTGCATAT
    sequence tr3: cc  =>   TGTGGTTAGCTGGTTGCATATCC
    sequence tr3: AT  =>   TGTGGTTAGCTGGTTGCAT

* s5: offset nucleotides at the begining of the annotated miRNAs::


    precursor         => agcctgtggttagctggttgcatatcc
    annotated miRNA   =>     TGTGGTTAGCTGGTTGCATAT
    s5                => AGCCTGTG


* s3:offset nucleotides at the ending of the annotated miRNAs::
 

    precursor         =>  cctgtggttagctggttgcatatccgc
    annotated miRNA   =>    TGTGGTTAGCTGGTTGCATAT
    s3                =>                     ATATCCGC


* type: mapped on precursor or miRNA sequences
* ambiguity: number of different detected precursors

Example::

    seq			miRNA		start	end	mism	tr5	tr3	add	s5	s3	DB amb
    TGGCTCAGTTCAGCAGGACC    hsa-mir-24-2    50      67      0       qCC     0       0       0       0       precursor 1
    ACTGCCCTAAGTGCTCCTTCTG  hsa-miR-18a*    47      68      0       0       0       tG      ATCTACTG        CTGGCA  miRNA 1
