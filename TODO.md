- [X] add length enrichment function
- [X] implement bayes function to decide which sequences go where
- [X] sequence in the cluster - text
- [X] create modules to save sqlite3 results
- [X] add doc about web servers with downstream analysis:
** http://compbio.uthsc.edu/miR2GO/home.php
- [X] create small stats: # of sequences, # in seqs.ma, # in alignment, # in clusters
- [X] counts of aligned reads
- [X] load matrix ma after alignment parsed
- [X] NOW html to connect to sqlite3 db to viz results

- [X] implement tRNA-scan
- [ ] NOW implement small rna prediction: CoRaL is an option for everything http://wanglab.pcbi.upenn.edu/coral/  or
- [X] NOW change to python alignment tool
- [ ] implement mirna prediction http://evryrna.ibisc.univ-evry.fr/miRBoost/index.html
- [ ] implement submire target prediction module http://research.nhgri.nih.gov/software/SubmiRine/user_guide.shtml
- [ ] Add to report: https://github.com/pkerpedjiev/forna/

# modules

* QC: 
  - [ ] figures 5' positions by # unique sequences, this by class: mirna, trna ...
* report:
  - [ ] add (R file) mirna path enrichment to report
  - [X] limit number of sequences in html, like 1000 top


# miraligner

inputs will be fastq or bam file. Will need to download the mirbase annotation (and hairpin.fa) if not given (need version). Species is required (for now).

- [X] create data set to test module
- [X] align with razer3 if option on
- [X] fix alignement error: get correct positions and guess additions
- [X] decide best hits
- [X] compare to mirbase ref: t5 t3 add mism
- [ ] ouput in SAM/BAM with custom flag for isomiRs
- [X] make it compatible with old miraligner output and isomiRs


# downstream

- [ ] isomirViz: html code to browse miraligner results for a project: preparation with miraligner module
- [ ] seqclusteR: will be an object that make easier to retrieve information like counts, rownames, annotations ... tidy tables, and figures

# documentation

- [X] add install subcommand for mirbase installation/genomes
- [X] add html report link to documentation
- [ ] figure showing general pipeline steps: seqcluster-helper/bcbio-nextgen
- [X] add documentation for seqclusterViz html reader

