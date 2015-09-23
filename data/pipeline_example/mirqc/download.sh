mkdir -p raw
cd raw

wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE49nnn/GSE49816/matrix/GSE49816-GPL10999_series_matrix.txt.gz

for F in `zcat GSE49816-GPL10999_series_matrix.txt.gz  | grep Sample_supplementary_file_2 | cut -f 2- | tr '\t' '\n' | sed 's/\"//g'` ; do 
    echo $F; 
    wget -nH -r --cut-dirs=10 $F; 
done 

wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE49nnn/GSE49816/matrix/GSE49816-GPL17561_series_matrix.txt.gz

for F in `zcat GSE49816-GPL17561_series_matrix.txt.gz  | grep Sample_supplementary_file_2 | cut -f 2- | tr '\t' '\n' | sed 's/\"//g'` ; do 
    echo $F; 
    wget -nH -r --cut-dirs=10 $F; 
done 

fastq-dump *sra
