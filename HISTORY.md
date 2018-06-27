- 1.2.4a*

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
