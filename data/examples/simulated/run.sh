set -e

# seqcluster simulator --fasta origin.fa --out sim
bowtie -f -a --best --strata ~/soft/bcbio/genomes/Hsapiens/hg19/bowtie/hg19 sim.fasta -S sim.sam
samtools view -Sbh sim.sam > sim.bam
samtools sort  sim.bam sim_sort

rm -rf res sim_sort_*
rm -rf res
seqcluster cluster -a sim_sort.bam -m sim.ma -o res -r ~/soft/bcbio/genomes/Hsapiens/hg19/seq/hg19.fa
# seqcluster predict -j res/seqcluster.json -o predict --reference  ~/orch/groups/bcbio/genomes/Hsapiens/hg19/seq/hg19.fa
rm -rf res_bayes
seqcluster cluster -a sim_sort.bam -m sim.ma -o res_bayes -r ~/soft/bcbio/genomes/Hsapiens/hg19/seq/hg19.fa --method bayes
# seqcluster report -j res/seqcluster.json -o report -r  ~/orch/groups/bcbio/genomes/Hsapiens/hg19/seq/hg19.fa
