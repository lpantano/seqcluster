# create hsapiends annotation for small rna analisis
# cel235

FINAL="srna-transcripts.gtf"

if [ -e $FINAL ] ; then

    rm $FINAL

fi


if [ ! -e cel.gff3 ] ; then 

    wget ftp://mirbase.org/pub/mirbase/21/genomes/cel.gff3

fi

    awk '$3=="miRNA"' cel.gff3 | sed 's/=/ /g' >> $FINAL 

if [ ! -e Caenorhabditis_elegans.WBcel235.81.gtf.gz ] ; then

    wget ftp://ftp.ensembl.org/pub/release-81/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.81.gtf.gz

fi

    zcat  wgEncodeGencodeBasicVM4.txt.gz | grep -v exon  >> $FINAL



