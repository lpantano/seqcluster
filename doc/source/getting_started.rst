.. _getting_started:


***************
Getting started
***************


clustering of small RNA sequences
-------- 

seqcluster generates a list of clusters of small RNA sequences, where they map on the genome, and the abundance in all the sample of the project


**REMOVE ADAPTER**

I am currently using ``cutadapt``:
::
    cutadapt --adapter=$ADAPTER --minimum-length=8 --untrimmed-output=sample1_notfound.fastq -o sample1_clean.fastq -m 17 --overlap=8 sample1.fastq 

**COLLAPSE READS***
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

You should use an aligner to map seqs.fa to your genome. A possibility is bowtie. 
From here, we need a file in BAM format for the next step.
VERY IMPORTANT: the BAM file should be sorted

::

    bowtie -a --best --strata -m 5000 -f INDEX seqs.fa -S | samtools view -Sbh /dev/stdin | samtools sort -o /dev/stdout temp > seqs.sort.bam


or 

::

    STAR --genomeDir $star_index_folder --readFilesIn res/seqs.fastq --outFilterMultimapNmax 5000 --outSAMattributes NH HI NM --outSAMtype BAM SortedByCoordinate


**CLUSTERING**

::

    seqcluster cluster -a res/Aligned.sortedByCoord.out.bam  -m res/seqs.ma -b  $ANN_FILE1,$ANN_FILE2  -o res/clluster


* `-a` is the SAM file generated after mapped with your tool, which input has been seqs.fa
* `-m` the previous seqs.fa
* `-b` annotation files in bed format (see below examples)
* `-g` annotation files in gtf format (see below examples)
* `-i` genome fasta file used in the mapping step (only needed if -s active)
* `-o` output folder
* `-d` create debug logging
* `-s` construction of putative precursor (NOT YEP IMPLEMENTED)

Example of a bed file for annotation (the fourth column should be the name of the feature): 

::

    chr1    157783  157886  snRNA   0       -
Example of a gtf file for annotation (the third column should be the name of the feature): 

::

    chr1    source  intergenic      1       11503   .       +       .       .....

**OUTPUTS**

* counts.tsv: count matrix that can be input of downstream analyses
* seqcluster.json: json file containing all information
* run.log: all messages at debug level
* trace.log: to keep trace of algorithm decision
