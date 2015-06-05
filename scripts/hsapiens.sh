# create hsapiends annotation for small rna analisis

FINAL="hsapiens.gtf"

if [ ! -e hsa.gff3 ] ; then 

    wget ftp://mirbase.org/pub/mirbase/20/genomes/hsa.gff3

fi

    awk '$3=="miRNA"' hsa.gff3 | sed 's/=/ /g' >> $FINAL 

if [ ! -e wgRna.txt.gz ] ; then

    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/wgRna.txt.gz

fi

    zgrep -v 'hsa-' wgRna.txt.gz | awk '{print $2"\t.\tncrna\t"$3"\t"$4"\t.\t"$7"\t.\tname "$5";"}' >> $FINAL

if [ ! -e tRNAs.txt.gz ] ; then

    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/tRNAs.txt.gz
    
fi

    zcat tRNAs.txt.gz | awk '{print $2"\t.\ttRNA\t"$3"\t"$4"\t.\t"$7"\t.\tname "$5";"}' >> $FINAL

if [ ! -e rmsk.txt.gz ] ; then

    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz
    
fi

    zcat rmsk.txt.gz | awk '{print $6"\t.\trepeat\t"$7+1"\t"$8+1"\t.\t"$10"\t.\tname "$12";"}' >> $FINAL


if [ ! -e refGene.txt.gz] ; then

    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
    
fi

    zcat refGene.txt.gz | awk '{print $3"\t.\tgene\t"$5"\t"$6"\t.\t"$4"\t.\tname "$13";"}' >> $FINAL


