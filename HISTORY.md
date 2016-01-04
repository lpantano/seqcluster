- development

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
