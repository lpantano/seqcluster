.. _getting_started:


***************
Getting started
***************

Best practices are implemented in a `python framework`_.

.. _python framework: https://github.com/lpantano/seqcluster-helper/blob/master/README.md

clustering of small RNA sequences
-------- 

seqcluster generates a list of clusters of small RNA sequences, where they map on the genome, and the abundance in all the sample of the project


**REMOVE ADAPTER**

I am currently using ``cutadapt``:
::
    cutadapt --adapter=$ADAPTER --minimum-length=8 --untrimmed-output=sample1_notfound.fastq -o sample1_clean.fastq -m 17 --overlap=8 sample1.fastq 

**COLLAPSE READS**

To reduce computational time, I recommend to collapse sequences, also it would help to apply filters based on abundances.
Like removing sequences that appear only once.

::

   seqcluster collapse -f sample1_clean.fastq -o collapse

Here I am only using sequences that had the adapter, meaning that for sure are small fragments.

**PREPARE SAMPLES**
::
    seqcluster prepare -c file_w_samples -o res --minl 17 --minc 2 --maxl 45

the file_w_samples should have the following format:

::

	lane1_sequence.txt_1_1_phred.fastq      cc1
	lane1_sequence.txt_2_1_phred.fastq      cc2
	lane2_sequence.txt_1_1_phred.fastq      cc3
	lane2_sequence.txt_2_1_phred.fastq      cc4

two columns file, where the first column is the name of the file with the small RNA sequences for each sample, and the second column in the name of the sample.

The fastq files should be like this:

::

    @seq_1_x11
    CCCCGTTCCCCCCTCCTCC
    +
    QUALITY_LINE
    @seq_2_x20
    TGCGCAGTGGCAGTATCGTAGCCAATG
    +
    QUALITY_LINE
    </pre>

Where _x[09]  indicate the abundance of that sequence, and the middle number is the index of the sequence.

This script will generate: seqs.fa and seqs.ma. 
* seqs.fa: have unique sequences and unique ids
* seqs.ma: is the abundance matrix of all unique sequences in all samples

**ALIGNMENT**

You should use an aligner to map seqs.fa to your genome. A possibility is bowtie or STAR. 
From here, we need a file in BAM format for the next step.
VERY IMPORTANT: the BAM file should be sorted

::

    bowtie -a --best --strata -m 5000 -f INDEX seqs.fa -S | samtools view -Sbh /dev/stdin | samtools sort -o /dev/stdout temp > seqs.sort.bam


or 

::

    STAR --genomeDir $star_index_folder --readFilesIn res/seqs.fastq --alignIntronMax 1  --outFilterMultimapNmax 1000 --outSAMattributes NH HI NM --outSAMtype BAM SortedByCoordinate


**CLUSTERING**

::

    seqcluster cluster -a res/Aligned.sortedByCoord.out.bam  -m res/seqs.ma -g $GTF_FILE  -o res/cluster --db example


* `-a` is the SAM file generated after mapped with your tool, which input has been seqs.fa
* `-m` the previous seqs.fa
* `-b` annotation files in bed format (see below examples) [deprecated]
* `-g` annotation files in gtf format (see below examples) [recommended]
* `-i` genome fasta file used in the mapping step (only needed if -s active)
* `-o` output folder
* `-d` create debug logging
* `-s` construction of putative precursor (NOT YET IMPLEMENTED)
* `--db` (optional) will create sqlite3 database with results that will be used to browse data with html web page (under development)

Example of a bed file for annotation (the fourth column should be the name of the feature): 

::

    chr1    157783  157886  snRNA   0       -
    
Example of a gtf file for annotation (the **third** column should be the name of the feature and
the value after `gene name` attribute is the specific annotation): 

:: 

    chr1    source  miRNA      1       11503   .       +       .       gene name 'mir-102' ;


**REPORT**

This will create html report using the following command assuming the output of `seqcluster cluster` is at `res`::

	seqcluster report -j res/seqcluster.json -o report -r $GENONE_FASTA_PATH

where `$GENOME_FASTA_PATH` is the path to the genome fasta file used in the alignment.

**OUTPUTS**

* counts.tsv: count matrix that can be input of downstream analyses
* size_counts.tsv: size distribution of the small RNA by annotation group
* seqcluster.json: json file containing all information
* log/run.log: all messages at debug level
* log/trace.log: to keep trace of algorithm decision
* report/html/index.html: table with all clusters and the annotation with sorting option
* report/html/[0-9]/maps.html: `summary`_ of the cluster with expression profile, annotation, and all sequences inside
* report/html/[0-9]/maps.fa: putative precursor

.. _summary: https://rawgit.com/lpantano/seqcluster/master/data/examples_report/html/1/maps.html
