_ ..exmaples

***************
Examples of small RNA analysis
***************

miRQC data
-------

**about**


`mirRQC project <http://www.nature.com/nmeth/journal/v11/n8/full/nmeth.3014.html>`_

samples overview:

>> Universal Human miRNA reference RNA (Agilent Technologies, #750700), human brain total RNA (Life Technologies, #AM6050), human liver total RNA (Life Technologies, #AM7960) and MS2-phage RNA (Roche, #10165948001) were diluted to a platform-specific concentration. RNA integrity and purity were evaluated using the Experion automated gel electrophoresis system (Bio-Rad) and Nanodrop spectrophotometer. All RNA samples were of high quality (miRQC A: RNA quality index (RQI, scale from 0 to 10) = 9.0; miRQC B: RQI = 8.7; human liver RNA: RQI = 9.2) and high purity (data not shown). RNA was isolated from serum prepared from three healthy donors using the miRNeasy mini kit (Qiagen) according to the manufacturer's instructions, and RNA samples were pooled. Informed consent was obtained from all donors (Ghent University Ethical Committee). Different kits for isolation of serum RNA are available; addressing their impact was outside the scope of this work. Synthetic miRNA templates for let-7a-5p, let-7b-5p, let-7c, let-7d-5p, miR-302a-3p, miR-302b-3p, miR-302c-3p, miR-302d-3p, miR-133a and miR-10a-5p were synthesized by Integrated DNA Technologies and 5′ phosphorylated. Synthetic let-7 and miR-302 miRNAs were spiked into MS2-phage RNA and total human liver RNA, respectively, at 5 × 106 copies/μg RNA. These samples do not contain endogenous miR-302 or let-7 miRNAs, which allowed unbiased analysis of cross-reactivity between the individual miR-302 and let-7 miRNAs measured by the platform and the different miR-302 and let-7 synthetic templates in a complex RNA background. Synthetic miRNA templates for miR-10a-5p, let-7a-5p, miR-302a-3p and miR-133a were spiked in human serum RNA at 6 × 103 copies per microliter of serum RNA or at 5-times higher, 2-times higher, 2-times lower and 5-times lower concentrations, respectively. All vendors received 10 μl of each serum RNA sample.

**commands**

Data was download from GEO web with this `script <https://github.com/lpantano/seqcluster/blob/master/data/pipeline_example/mirqc/download.sh>`_. The following `config <https://github.com/lpantano/seqcluster/blob/master/data/pipeline_example/mirqc/config.txt>`_  file was used with `seqcluster-helper <http://seqcluster.readthedocs.org>`_ to analysis the data with the following command::
  
  seqcluster-helper.py --sample-map config.txt  --db ~/groups/seqcluster/data/mirbase --aligner-index ~/groups/bcbio/genomes/Hsapiens/hg19/star --gtf-file ~/groups/seqcluster/data/annotation/hsapiens.gtf --species hsa --adapter TGGAATTCTCGGGTGC --reference ~/groups/bcbio/genomes/Hsapiens/hg19/seq/hg19.fa`.

* Note: `bcbio <http://github.com/chapmanb/bcbio-nextgen>`_ will wrap seqcluster-helper in the next month. Follow @seqbuster to know the exact date.
