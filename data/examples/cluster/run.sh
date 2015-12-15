seqcluster cluster -a seqs_map.bam -m seqs_set.ma -g ann_reduced.gtf -o test_out_cluster --db example -ref ../../genomes/genome.fa

seqcluster report -j test_out_cluster/seqcluster.json -o test_out_report -r ../../genomes/genome.fa

# rm -rf predict; seqcluster predict -j res/seqcluster.json --bed res/positions.bed -o test_predict -r ../genomes/genome.fa --coral --bam seqs_map_rmlw.bam
