X add length enrichment function
X implement bayes function to decide which sequences go where
X sequence in the cluster - text
* add doc about web servers with downstream analysis:
** http://compbio.uthsc.edu/miR2GO/home.php
* create small stats: # of sequences, # in seqs.ma, # in alignment, # in clusters
* segmentation algorithm: RSEQtools has a C tool to create segmentation: https://github.com/sbonerlab/rseqtools/blob/master/src/bgrSegmenter.c
* NOW implement tRNA-scan
* Add to report: https://github.com/pkerpedjiev/forna/
* implement small rna  prediction: CoRaL is an option for everything http://wanglab.pcbi.upenn.edu/coral/  or
* implement mirna prediction http://evryrna.ibisc.univ-evry.fr/miRBoost/index.html
* implement submire target prediction module http://research.nhgri.nih.gov/software/SubmiRine/user_guide.shtml
* counts of sequences - barplot / problem: need to be in the same page, select and update :( clickId?? maybe
* annotation - datatable
* viz with http://stanford.edu/~mwaskom/software/seaborn/index.html or D3.js
* trace: tree of initial cluster and how ends up to a final cluster/s. To plot in R needs threen columns file: 
** loci1 loci2 0 %similarity => merge cluster
** loci1 newloci nseqs 0 => iter_loci
** loci2 newloci nseqs 0 => iter_loci
