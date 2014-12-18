After the comparison of tools to annotate miRNAs explained here, I decided to try STAR ([which was one of the best](http://lorenapantano.wordpress.com/2014/02/28/mirna-annotation-tools-which-is-the-best/)) to map miRNA sequences against the whole genome instead against miRBase to reduce the number of step in the pipeline.

The step in the pipeline are:

* simulate miRNA sequences with SeqBuster simulator
* index human genome
* map with STAR
* [optional] filter not primary hits with bamtools
* intersect with bedtools using hsa.gff3 file
* parse results with python script

The [result](https://raw.githubusercontent.com/lpantano/reproducibility/master/mirannotation/star_genome/stats_star.html) show how is not good to filter out secondary hits, and the accuracy is reduced dramatically using whole genome for mapping.

![annotation of different isomiRs using STAR](https://raw.githubusercontent.com/lpantano/reproducibility/master/mirannotation/star_genome/stats_star_files/figure-html/iso-no-filtered-1.png)

It misses mainly sequences with mismatches and additions at the 3' end.

Highly recommend to map against miRNA first in this case, and then against the genome to decide best hit.
