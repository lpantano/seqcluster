After the comparison of tools to annotate miRNAs explained here, I decided to try [STAR](http://code.google.com/p/rna-star/) ([which was one of the best](http://lorenapantano.wordpress.com/2014/02/28/mirna-annotation-tools-which-is-the-best/) using miRBase database) to map miRNA sequences against the whole genome instead against [miRBase](http://www.mirbase.org/) to reduce the number of steps in the pipeline.

The steps in the pipeline are:

* simulate miRNA sequences with [SeqBuster](http://github.com/lpantano/seqbuster) simulator
* index human genome
* map with STAR
* [optional] filter not primary hits with bamtools
* intersect with bedtools using hsa.gff3 file
* parse results with python script

The [result](https://rawgit.com/lpantano/reproducibility/master/mirannotation/star_genome/stats_star.html) shows how is not good to filter out secondary alignments, and furthermore, the accuracy is reduced dramatically using whole genome directly when mapping. In that link you will get the commands used for the analysis.

![annotation of different isomiRs using STAR](https://raw.githubusercontent.com/lpantano/reproducibility/master/mirannotation/star_genome/stats_star_files/figure-html/iso-no-filtered-1.png)

From the 16900 initial sequences, 16377 mapped to the genomea, and only 12549 overlapped with any miRNA. The correctly annotated sequences are the ones without mismatches or additions, meaning that you will lose miRNA sequences with post-transcriptional modifications.

Highly recommend to map against miRNA first in this case, and then against the genome to decide best hit.
