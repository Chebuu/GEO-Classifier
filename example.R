###############################
### EXAMPLE IMPLEMENTATION  ###
###############################
source('methods.R')
RUN_LOCAL <- F # Has a local dataset been saved that should be used instead of querying GEO?
SAVE_DATASET <- F # Write the dataset to a file?
RETRY_ON_ERROR <- T # If a GEO request fails, should the script indefinitely continue trying to download the files?
LOCAL_DATA <- 'data.txt' # Where is local data saved? Or, where should it be saved after querying GEO? This is a relative path.

# Choose some GEO datasets (GSEs) and inlclude class labels for each sample (GSM) in the dataset.
gse.list <- list('GSE72078', 'GSE22792', 'GSE35028', 'GSE35029', 'GSE32658', 'GSE28633', 'GSE39876')
labels.list <- list(c(rep(3, length(1:8)), rep(1, length(9:12))), # ----------------------------------------------------------------------------------- # GSE72078
                    c(rep(3, length(1:8)), rep(1, length(9:10)),rep(2, length(11:14))), # ------------------------------------------------------------- # GSE22792
                    c(rep(1, length(1:6)), rep(2,length(7)), rep(3,length(8:17))), # ------------------------------------------------------------------ # GSE35028
                    c(rep(2, length(1:3)), rep(1, length(4:30)), rep(1, length(31:32)),rep(0,length(33:36)),rep(2,length(37:37)),rep(3,length(38:47))), # GSE35029
                    c(rep(2, length(1:72))), # -------------------------------------------------------------------------------------------------------- # GSE32658
                    c(rep(2, length(1:3)), rep(0, length(4:6)), rep(1, length(7:12))), # -------------------------------------------------------------- # GSE28633
                    c(rep(2, length(1:6)), rep(1, length(7:12)))) # ----------------------------------------------------------------------------------- # GSE39876

# The samples to withhold from training and use for testing
# Redundant nested vectors here are simply used to demarcate the true class of each sample
class.charlevels <- c('Adult', 'ESC', 'iPSC') # [1=Adult, 2=ESC, 3=iPSC] Character labels are needed by caret::train to compute class probabilities
test.cols <- c(c('GSM860992', # Adult keratinocytes (GSE35029)
                 'GSM860995', # Peripheral CD34+ progenitors (GSE35028)
                 'GSM980629', # Neural progenitors (GSE39876)
                 'GSM563510'),# MRC5 fibroblasts (GSE22792)
               c('GSM860962', # Embryonic stem cells (GSE35027)
                 'GSM980627', # Embryonic stem cells (GSE39876)
                 'GSM980628', # Embryonic stem cells (GSE39876)
                 'GSM709608'),# Embryonic stem cells (GSE28633)
               c('GSM860998', # iPSCs derived from cord blood (non-viral) (GSE35029)
                 'GSM861006', # iPSCs derived from fibroblasts (viral) (GSE35029)
                 'GSM563502'))# iPSCs derived from culture with Yaminaka factors (GSE22792)


if(!RUN_LOCAL){
  dataset.labeled <- makeLabeledGEODataset(gse.list,labels.list, retry.on.error=RETRY_ON_ERROR)
  if(SAVE_DATASET){
    write.table(dataset.labeled, file=LOCAL_DATA, row.names=T, col.names=T, append=F)
  }
}else{
  dataset.labeled <- read.table(LOCAL_DATA, header=T, row.names = 1)
}
features.subset <- subsetFeatures(dataset.labeled, n.features=1000, plot.MDS=T)
res <- trainAndTestClassifier(dataset.labeled, features.subset=features.subset, class.charlevels=class.charlevels, test.cols=test.cols)
print(res)