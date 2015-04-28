"""Use CoRaL to predict small RNA function"""

### Running CoRAL
## call intervals corresponding to discrete small RNA-producing loci
# produces coral/loci.bed
# call_smrna_loci.sh $bam $conf

## compute features
# feature_lengths.sh $bam $conf
# feature_antisense.sh $bam $conf
# feature_entropy.sh $bam  $conf $annot/chromInfo.txt
# feature_nuc.sh $bam $conf
# feature_mfe.sh $bam $conf $annot/hsa19.fa $annot/chromInfo.txt

## label the loci based on known annotation data - this is only needed for training
# annotate_loci.sh coral/loci.bed  $annot/hsa19.gff $annot/class_pri.txt

## generate data_x.txt and data_y.txt for input into the training and/or prediction
# make_data_matrix.sh coral

## train a random forest classifier on the generated features for 3 classes
# the result will be in coral/run_xxxxx where xxx is a hash on the parameters used
# parameters used here: require 15 reads at a locus, use these three classes only
# coral_train.R -r 15 -c "miRNA,snoRNA_CD,tRNA" \
#          coral/data_x.txt coral/data_y.txt coral

## predict on entire dataset and use known data (data_y) to assess training performance
## and the model that was trained and outputted to "coral/run_*/"
# places results in pred_out dir
# coral_predict.R -r 15 -y coral/data_y.txt coral/data_x.txt \
#         coral/run_* pred_out

def cmd(subcmd):
    cmd = ""
    return cmd
