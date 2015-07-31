# create hsapiends annotation for small rna analisis
# mm10

FINAL="mmusculus.gtf"

if [ -e $FINAL ] ; then

    rm $FINAL

fi


if [ ! -e mmu.gff3 ] ; then 

    wget ftp://mirbase.org/pub/mirbase/20/genomes/mmu.gff3

fi

    awk '$3=="miRNA"' mmu.gff3 | sed 's/=/ /g' >> $FINAL 

if [ ! -e wgEncodeGencodeBasicVM4.txt.gz ] ; then

    wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/wgEncodeGencodeBasicVM4.txt.gz

fi

    zcat  wgEncodeGencodeBasicVM4.txt.gz | awk '{print $3"\t.\tencode\t"$5"\t"$6"\t.\t"$4"\t.\tname "$13";"}' | awk '$5-$4 < 500' >> $FINAL

if [ ! -e tRNAs.txt.gz ] ; then

    wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/tRNAs.txt.gz
    
fi

    zcat tRNAs.txt.gz | awk '{print $2"\t.\ttRNA\t"$3"\t"$4"\t.\t"$7"\t.\tname "$5";"}' >> $FINAL

if [ ! -e rmsk.txt.gz ] ; then

    wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/rmsk.txt.gz
    
fi

    zcat rmsk.txt.gz | awk '{print $6"\t.\trepeat\t"$7+1"\t"$8+1"\t.\t"$10"\t.\tname "$12";"}' >> $FINAL


if [ ! -e refGene.txt.gz ] ; then

    wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/refGene.txt.gz
    
fi

    zcat refGene.txt.gz | awk '{print $3"\t.\tgene\t"$5"\t"$6"\t.\t"$4"\t.\tname "$13";"}' >> $FINAL


