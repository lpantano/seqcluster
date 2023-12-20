- 1.2.10

  * Fix bug counting sequences expression in stats file.

- 1.2.9

  * Fix bug where UMI is mistakenly detected in read names containing "ILLUMINA"

- 1.2.8

  * Fix bug when detecting the end of the cluster due to big gaps in biopython function

- 1.2.7

  * Fix bug when writing files for debug of big meta-clusters: 
    https://github.com/lpantano/seqcluster/issues/47
    https://github.com/bcbio/bcbio-nextgen/issues/2948
  * Add version option
  
- 1.2.5

  * Fix error when the precursor is too long to ignore RNAfold calculation. Thanks to @ZhuZhuoHSPH and @kthlnktng

- 1.2.4
  
  * Fix multiple errors when running in python 3 due to map function.
  * Fix error in collapsing fasta files.
  * Fix end of line character for counts_sequence.tsv.
  * Remove `map` function from quality class in collapse function to avoid seg.fault in python3.
  * Use DESeq2 normalization strategy.
  * Fix more errors in python3 env.
  * Fix UMI checking when the input file is a gzip file. Thanks to @rbatorsky-claritas.
  * Fix header bug
  * Initiate migration to py3*. Thanks to @smoe.
  * Include mirtop annotation.
  * Fix upgrade command.
  * Fix UMIs detection to count using unique seq + umi. Thanks to @mshadbolt
  * Remove Cpy code and use biopython
  * Clean test examples
  * Fix UMI error when sequences have different sizes. Thanks to @mshadbolt
  * Support UMI tag when collapsing
  * Add count matrix for each sequence
  * Remove HTML report
  * Allow size parameter during collapsing reads
  * Fix reporting DB when precursor is masked
  * Add conflict to output
  * Fix bug in prepare sample that will setup min-shared
    to samples size always.

- 1.2.3

  * Add --feature_id as an option to specify the attribute
    to use in the GTF file for annotation
  * Add gene_id as a 2nd option to add GTF annotation
  * Only do rnafold for precursors shorter than 200nt

- 1.2.2

  * Use bedtools for bamtobed and clustering
  * Only update seqcluster code when upgrading
  * Add beter logging to prepare sub-command

- 1.2.1

  * Fix expression profile when no sequence at that position
  * Fix reading from profile file to avoid calculation

- 1.2.0

  * Add function to map SNPs to genome coordinate
  * Add RNAfold to report for html vis.
  * Adapt C code to macosx
  * Improve test functions

- 1.1.14

  * Improve miRNA annotation function and add first functions
  to allow SNP detections
  * remove bcbio funtions to simpler to avoid circular dependency
  * add function to get targets from targetscan for human
  * select best precursor based on size
  * parse pyMAtch to be compatible with other tools
  * add C function working as miraligner by @franpantano (pyMatch)

- 1.1.11

  * first version of simulator of 33 nt long sRNA from fasta
  * add Rmd template for bcbio pipeline
  * add peak detection to cluster, as a try to get the mature sequence

- 1.1.8

  * fix indexing bam file error when running under python
  * remove colorlog dependency
  * create light package in conda for bcbio installation
