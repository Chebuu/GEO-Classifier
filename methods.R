####################
##  Description   ##
####################
# This script contains a collection of methods for the classification of RNAseq data using the caret package

##############
##   TODO   ##
##############

# Unecessary class.charlevels
## There is no need for the class.charlevels parameter
## The user can just rep('myClass', length(1:10)) or use integers which can be coerced to strings

# Support different GPLs
## Right now, using different GPLs will result in 0 or few genes after subsetting
## The fix would be to convert all probes to ENTREZ IDs and subset using those IDs

# Metadata column selection of class labels
## If the item in labels.list is a string, assume this is the metadata column in the GSE that contains class labels
## Will need to icnlude a way to match these labels to other labels passed by the user.

## [DONE] Support named gse.list arguments
# gse.list and labels.list (args to makeLabeledGEODataset) should be combined into a single named list.
# Each vector (element) in the list can be either numerical or a single string specifying the metadata column to use for labels
## makeLabeledGEODataset(gse.list=list(GSE72078=c(rep(3, length(1:8)), rep(1, length(9:12))), labels.list=NULL, ...)
## if gse.list is not named, labels.list should not be null, and if it is throw an error, because we can't find the datasets without GSE identifiers

library(rlang)
library(GEOquery)
library(edgeR)
library(caret)
library(colorspace)

source('GEOfix.R')

###############
### METHODS ###
###############

retryOnError <- function(func, args, n.tries=Inf) {
  # TODO: The call to func(args) cannot return FALSE, because the program will spin indefinitely. This is a problem.
  # Repeatedly tries to call a function until the call does not throw an error.
  # Returns the result of do.call(func, args)
  # @ param {func} A function definition to be called with args.
  # @ param {args} A list of arguments to be passed to func. May be a named list. 
  # @ param {n.tries} An integer specifying the number of times do.call(func, args) should be allowed to throw an error before giving up and returning FALSE.
  # @ example retryOnError(do.call, list(c,1)) # Executing do.call(c,1) will throw an error, because the second argument to do.call must be a list.
  # @ example retryOnError(makeLabeledGEODataset, list(gse.list, labels.list))
  tempenv <- new.env()
  tempenv$triesLeft <- n.tries
  errorCallBack <- function(e){
    print(e); print('Retrying...') 
    tempenv$triesLeft <- tempenv$triesLeft - 1
    return(F)
  }
  repeat{
    res <- tryCatch(do.call(func, args), error=errorCallBack)
    if(!is_false(res) | tempenv$triesLeft == 0){
      if(!is_false(res)) {
        print('Success.')
      }else {
        print(sprintf('Maximum tries (%s) reached.', n.tries))
      }
      break
    }
  }
  return(res)
}

argmax <- function(Y) {
  # Returns the column index of the maximum value in each row.
  # @ param {Y} A numeric matrix of model outputs. 
  apply(Y, 1, function(y){
    which(y == max(y))
  })
}

makeLabeledGEODataset <- function(gse.list, labels.list=NULL, n.tries=Inf, normalize.counts=edgeR::cpm, ...) {
  # Returns a labeled dataset that can be used to train a classifier. The labels exist in the first row of the returned dataset. The features of this dataset are the common subset of features (probes) among each individual expression set. Note, the dimensions of this dataset are probesXsamples, which is the transpose of the canonical machine learning dataset (samplesXfeatures).
  # @ param {geo.list} A list of GSE accession numbers as strings. All GSEs should be of the same GPL to ensure that a common subset of genes exists among each GSE in the list. To include individual GSMs, simply include the GSM accession number and a single label for the GSM in labels.list.
  # @ param {labels.list} A list of integer labels for the samples (columns) of each gse.list member. Each labels.list member corresponds to one GSE or GSM in gse.list. In other words, labels.list holds the labels that will be appended to the dataset and used for classification. To exclude a particular sample from a GSE, label it as 0 and it will be eliminated. It is important to note that the length of labels.list must be equal to the total number of GSMs (smaples) acquired from a query to GEO using gse.list, and that the order of labels in labels.list must match the order of GSMs (samples) in gse.list. 
  # @ param {n.tries} The number of times to retry GEO queries. Sometimes a query will error with 404 status, and simply resending the query will succeed with 200 status. 
  # @ param {normalize.counts} The normalization function to apply to each individual GSE in gse.list. The default is counts-per-million. Pass NULL argument to skip normalization and use raw counts.
  # @ param {...} Additional arguments to be passed to caret::train()
  # @ example makeLabeledGEODataset(gse.list=list('GSE35029'), labels.list=list(c(rep(2, length(1:3)), rep(1, length(4:30)))))
  # @ example makeLabeledGEODataset(list('GSE35029'= list(c(rep(2, length(1:3)), rep(1, length(4:30)))))
  
  reduce.genes <- function(GSE){
    if(length(GSE@header$platform_id) > 1){
      # If there is more than one platform in this GSE, then there will be GSMs with different probes
      print(sprintf('Subsetting genes in %s', GSE@header$geo_accession))
      genes.common <- Reduce(intersect, lapply(GSE@gsms, function(gsm){ # Find the common subset of genes among all GSMs in the GSE
        # TODO: should just Reduce(intersect over each GPL, not each GSM (but this block is relatively fast as is)
        return(gsm@dataTable@table[,1]) # Get the genes from each GSM in the GSE
      }))
      GSE@gsms <- lapply(GSE@gsms, function(gsm){ # Loop over all GSMs in the GSE
        print(sprintf('%s', gsm@header$geo_accession))
        genes.common.idx <- which(gsm@dataTable@table[,1] %in% genes.common) # Doing this by index speeds things up considerably
        gsm@dataTable@table <- gsm@dataTable@table[genes.common.idx,] # Return only the common subset of genes from each GSM
        return(gsm)
      })
      print('Complete')
    }
    return(GSE)
  }
  
  GSM2matrix <- function(GSM) {
    # Lambda function used by GSE2matrix to convert a given GSM to a 1 X N matrix, with rownames as gene probes
    m <- matrix(GSM@dataTable@table[,2]) # Dim 2 is the count vector
    rownames(m) <- GSM@dataTable@table[,1] # Dim 1 is the vector of probe names
    return(m)
  }
  
  GSE2matrix <- function(GSE) {
    # Lambda function used to convert the GSEs in the dataset to matrices
    GSE <- reduce.genes(GSE) # Some GSMs under a single GSE have different number of rows. Return only the common rows in order to form a rectangular matrix.
    m <- do.call(cbind, lapply(GSE@gsms, GSM2matrix)) # Create the matrix
    colnames(m) <- names(GSE@gsms) # Column names are sample names in this case (not features)
    return(m)
  }
  
  if(is.null(n.tries))
    n.tries <- FALSE
  
  if(is.null(normalize.counts) | !isTRUE(normalize.counts))
    normalize.counts <- function(x) return(x) # If the user passes normalize.counts=FALSE, do not normalize the count matrices
  if(isTRUE(normalize.counts) & !is.function(normalize.counts))
    normalize.counts <- edgeR::cpm # If the user errantly passes normalize.counts=TRUE, just reset the normalization function to the default
  
  # This block adds support for using a named gse.list instead of separating the two into gse.list and labels.list
  if(length(names(gse.list)) > 0){
    # If gse.list is a named list, then...
    if(length(names(gse.list)) == length(gse.list)){
      # If all items in the list are named, then ...
      labels.list <- gse.list # Get the class labels 
      gse.list <- names(gse.list) # Get the GSE identifiers
    }else{
      stop('Please name all list items in gse.list')
    }
  }else if(length(labels.list) == length(gse.list)){
    # The user did not name gse.list, rather they included labels.list, which is the default expectation
    # This could be refactored as length(labels.list) != length(gse.list) and then throw the error, but I think this is more clear.
    NULL # Break
  }else{
    stop('Unable to infer labels for gse.list. Please check args to gse.list and labels.list .')
  }
  
  gse.all <- lapply(gse.list, function(GSE){ 
    # Query GEO with retries if the query throws an error (e.g. 404 status)
    if(is_false(n.tries)){
      getGEO.simple(GSE)
    }else{
      retryOnError(getGEO.simple, list(GSE), n.tries)  
    }
  })
  
  # TODO: The variable name gse.all.subset can be changed to gse.all to overwrite gse.all, because gse.all is not used again after gse.all.subset is declared. In other words, there is no reason to declare two variables.
  gse.all <- lapply(gse.all, GSE2matrix) # Convert all @dataTables into matrices so normalization and other downstream methods work as expected
  gse.all <- lapply(gse.all, function(gse){ normalize.counts(gse, ...) }) # Apply the normalization function to each individual GSE
  genes.common <- Reduce(intersect, lapply(gse.all, rownames)) # Find the common subset of genes among all GSMs from all GSEs in the dataset
  gse.all.subset <- lapply(gse.all, function(x){ x[genes.common,]}) # Inclue only the common subset of genes
  gse.all.subset <- do.call(cbind, gse.all.subset) # Combine all the data sub-sets into a single matrix 
  gse.all.subset.labeled <- rbind(unlist(labels.list), gse.all.subset) # Append class labels
  if(any(unlist(labels.list) == 0)) gse.all.subset.labeled <- gse.all.subset.labeled[,-(which(unlist(labels.list) == 0))] # Remove GSMs with a class label of 0
  return(gse.all.subset.labeled)
}

subsetFeatures <- function(dataset.labeled, n.features=NULL, plot.MDS=T) {
  # Returns the top n.features of differentially expressed genes among the classes in dataset.labeled
  # @ param {dataset.labeled} The complete dataset created with a call to makeLabeledGeoDataset(), the first row of which is a vector of class labels.
  # @ param {n.features} An integer describing the desired number of features to keep. This will select from a list of differentially expressed genes sorted by p-value. Default is to select the top 50% of differentially expressed genes.
  # @ param {plot.MDS} Include a multidimentional scaling plot of log2 fold changes with a call to limma::plotMDS()
  # @ example subsetFeatures(makeLabeledGEODataset(list('GSE35029'), list(c(rep(2, length(1:3)), rep(1, length(4:30))))))
  if(is.null(n.features)) n.features <- floor((nrow(dataset.labeled)-1)/2)
  classes <- dataset.labeled[1,]
  data <- as.matrix(dataset.labeled[-(1),], ncol=ncol(dataset.labeled))
  cols <- rainbow(nlevels(as.factor(unlist(classes))))[unlist(classes)]
  model <- model.matrix(~unlist(classes))
  dlist <- DGEList(counts=data, group=unlist(classes), remove.zeros = T)
  dlist <- estimateDisp(dlist, design=model)
  fit <- eBayes(lmFit(dlist$counts, design=model))
  features.subset <- rownames(topTable(fit, number=n.features))
  if(plot.MDS) plotMDS(dlist, labels=classes, col=cols)
  return(features.subset)
}

trainAndTestClassifier <- function(dataset.labeled, features.subset, test.cols=NULL, class.charlevels=NULL, method.train='rbfDDA', method.cv='adaptive_cv', caret.trainControl=NULL, ...) {
  # TODO: Include elipses for further arguments to be passes to caret::train . For example, preprocessing. Maybe just make a custom caret method.
  # TODO: Change method.train and method.cv to just model and cv
  # Trains a caret classifier with labeled training data
  # Returns a list
  # @ param {dataset.labeled} A labeled dataset containing both training and test data. The first row of which is a vector of class labels. The rows of the dataset are features, and the columns of the dataset are samples. 
  # @ param {features.subset} A vector of rownames (not column names) or indices representing the subset of features in dataset.labeled to be used to train the classifier. If features.subset is a vector of row indices, those indices should be relative to the dataset after the training labels are removed from the first row of dataset.labeled. A character vector of rownames is recommended for this reason. It is important to remember that dataset.labeled is a matrix of gene expression values, where features (genes) are represented as rows, and samples (cells/runs) are represented as columns. This may be counter intuitive, because data is typically formatted as the transpose, with columns representing features and rows representing samples.
  # @ param {class.charlevels} A character vector whose items denote the human-readable meaning of each numeric class label. For example, if dataset.labeled[1,] == c(3,3,2,1,2,3,) then class.charlevels might be c('one','two','three') and an integer label of 1 will be correctly given the character label of 'one'. Another example, for the same class labels, class.charlevels could be c('Adult', 'ESC', 'iPSC').
  # @ param {test.cols} A character vector of column names that should be removed from the training dataset and used to test the model. Integer column indices are also accepted.
  # @ param {method.train} A string specifying any of the acceptable caret classification/methods methods. It will simply be passed as the "method" argument to caret::train(). From the caret::train() man page: A string specifying which classification or regression model to use. Possible values are found using names(getModelInfo()). See http://topepo.github.io/caret/train-models-by-tag.html. A list of functions can also be passed for a custom model function. See http://topepo.github.io/caret/using-your-own-model-in-train.html for details.
  # @ param {method.cv} A string specifying any of the acceptable caret cross-validation methods. It will simply be passed as the "method" argument to caret::trainControl(). From the caret::trainControl() man page: The resampling method: "boot", "boot632", "optimism_boot", "boot_all", "cv", "repeatedcv", "LOOCV", "LGOCV" (for repeated training/test splits), "none" (only fits one model to the entire training set), "oob" (only for random forest, bagged trees, bagged earth, bagged flexible discriminant analysis, or conditional tree forest models), timeslice, "adaptive_cv", "adaptive_boot" or "adaptive_LGOCV"
  
  classes <- unlist(dataset.labeled[1,]) # The first row of dataset.labeled contains the class labels
  dataset <- dataset.labeled[-(1),] # subsetFeatures returns row names, but labels are removed from the dataset incase keep.features is passed as a vector of integers.
  dataset <- dataset[features.subset, ] # Can't mix negative subscripts and indices.
  
  if(is.null(test.cols)) test.cols <- sample(1:ncol(dataset), round(ncol(dataset)*0.1))
  if(is.character(test.cols)) test.cols <- which(colnames(dataset) %in% test.cols)
  if(!is.character(classes)) classes.char <- class.charlevels[classes]
  
  x.train <- dataset[ ,-(test.cols)]
  y.train.int <- classes[-(test.cols)]
  y.train.char <- classes.char[-(test.cols)]
  x.test <- dataset[ ,test.cols]
  y.test <- classes[test.cols]
  y.test.char <- classes.char[test.cols]
  
  which.duplicate <- 'No test samples are present in the training set'
  if(any(colnames(x.test) %in% colnames(x.train))){
    which.duplicate <- colnames(x.test) %in% colnames(x.train)
    print(which.duplicate)
    warning('The above GSMs are present in both the training and test sets. This may cause overfitting.')
  }
  
  if(is.null(caret.trainControl)){
    seeds <- vector(mode = 'list', length = 51)
    for(i in 1:50) seeds[[i]] <- sample.int(1000, 27)
    seeds[[51]] <- sample.int(1000, 1)
    ctrl <- trainControl(method=method.cv, classProbs=T, seeds=seeds)  
  }else{
    ctrl <- caret.trainControl
  }
  
  model <- caret::train(x=t(x.train), y=as.factor(y.train.char), method=method.train, trControl=ctrl, ...)
  model.argmax <- argmax(model$finalModel$fitted.values)
  model.confMat <- caret::confusionMatrix(
    as.factor(model.argmax),
    as.factor(unlist(y.train.int))
  )
  
  predictions <- predict(model, t(x.test))
  predictions.int <- match(predictions, class.charlevels)
  predictions.confmat <- caret::confusionMatrix(
    as.factor(predictions), 
    as.factor(unlist(y.test.char))
  )
  
  return(
    list(
      model.train=list(
        model=model,
        fitted.values=model.argmax,
        predictions.int=y.train.int,
        confusion.matrix=model.confMat
      ),
      model.test=list(
        predictions.char=predictions,
        predictions.int=predictions.int,
        actual.char=unlist(y.test.char), 
        actual.int=unlist(y.test),
        confusion.matrix=predictions.confmat
      ),
      misc=list(
        which.duplicate=which.duplicate
      )
    )
  )
}
