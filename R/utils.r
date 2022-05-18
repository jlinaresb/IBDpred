# Utils

# This function select those samples from the same patients and select one by last timestamp!!
# Arguments:
#   - data: clinical data
# Return:
#   - resData: clinical data with length(unique(resData$HOST_SUBJECT_ID)) == dim(resData)[1]
selectDupPatients = function(data){
  
  # Check
  stopifnot('HOST_SUBJECT_ID' %in% names(data))
  stopifnot('COLLECTION_TIMESTAMP' %in% names(data))
  
  dd = duplicated(data$HOST_SUBJECT_ID)
  ids = data$HOST_SUBJECT_ID
  duplicates = unique(ids[which(dd == T)])
  
  toDelete = list()
  for (i in 1:length(duplicates)) {
    
    pats = data[which(data$HOST_SUBJECT_ID == duplicates[i]), ]
    maxDate = rownames(pats[which(pats$COLLECTION_TIMESTAMP == max(pats$COLLECTION_TIMESTAMP)), ])
    if (length(maxDate) == dim(pats)[1]) {
      maxDate = maxDate[1]
    }
    toDelete[[i]] = setdiff(rownames(pats), maxDate)
  }
  del = unlist(toDelete)
  toStay = setdiff(rownames(data), del)
  res = data[match(toStay, rownames(data)), ]
  
  # Remove those patients still duplicated (same sample, so I retrieve randomless!)
  
  DUP = res$HOST_SUBJECT_ID[which(duplicated(res$HOST_SUBJECT_ID) == T)]
  random = list()
  for (i in 1:length(DUP)) {
    kk = res[which(res$HOST_SUBJECT_ID == DUP[i]), ]
    select = rownames(kk)[1]
    
    random[[i]] = setdiff(rownames(kk), select)
  }
  rm = unlist(random)
  sel = setdiff(rownames(res), unlist(random))
  dataRes = res[match(sel, rownames(res)), ]
  
  return(dataRes)
}


splitData = function(dirData, file, outDir, save = TRUE, keep = FALSE){

  phyloseq = readRDS(paste0(dirData, file))

  # Access to otu_table
  otu = t(otu_table(phyloseq)@.Data)
  otu = as.data.frame(otu)
  rownames(otu) = make.names(rownames(otu))

  # Delete those variables with all zero values
  zeros = apply(otu, 2, function(x) sum(x))
  noZeros = which(zeros != 0)
  otu = otu[, noZeros]

  # Save in a vector the colnames and rownames of otu table
  cols = tax_table(phyloseq)@.Data[,6]
  i = match(colnames(otu), names(cols))
  cols = cols[i]
  rows = make.names(rownames(otu))

  # Y variable
  sampleData = sample_data(phyloseq)

  Yvars = c("ages", "AGE_YEARS", "SUBSET_ANTIBIOTIC_HISTORY")
  Yvars = sampleData[, intersect(Yvars, names(sampleData))]

  # Covariates
  cvrts = c("COUNTRY",                                                                        
            "ALCOHOL_CONSUMPTION",  
            "EXERCISE_FREQUENCY",                                                                      
            "DIET_TYPE",                                                                     
            "ALCOHOL_FREQUENCY",                                                                       
            "SMOKING_FREQUENCY",                                                                       
            "SLEEP_DURATION",                                                                          
            "ANTIBIOTIC_HISTORY",                                                                      
            "SUBSET_ANTIBIOTIC_HISTORY",                                                               
            "HEIGHT_CM",                                                                               
            "BMI",                                                                                     
            "WEIGHT_KG",                                                                               
            "SUBSET_HEALTHY",                                                                          
            "APPENDIX_REMOVED",
            "RACE_African American",
            "RACE_Asian or Pacific Islander",
            "RACE_Caucasian",
            "RACE_Hispanic",
            "SEX_male",
            "LIVER_DISEASE_yes",
            "IBD_yes",
            "IBS_yes",
            "DIABETES_yes",
            "CDIFF_yes",
            "CARDIOVASCULAR_DISEASE_yes",
            "LUNG_DISEASE_yes",
            "AUTOIMMUNE_yes")


  cvrts = make.names(cvrts)

  ClinCvrts = sampleData[, intersect(cvrts, names(sampleData))]
  #otu = cbind.data.frame(cvrts, otu)
  rownames(ClinCvrts) = make.names(rownames(ClinCvrts))
  stopifnot(rownames(otu) == rownames(ClinCvrts))

  # Split data into train y test (80% to train)
  set.seed(2)
  split = sample(nrow(otu), floor(0.8*nrow(otu)))
  OtuTrain = otu[split,]
  OtuTtest = otu[-split,]

  ClinTrain = ClinCvrts[split,]
  ClinTest = ClinCvrts[-split,]

  YTrain = Yvars[split,]
  YTest = Yvars[-split,]


  res = list(OtuTrain = OtuTrain, OtuTest = OtuTtest,
             ClinTrain = ClinTrain, ClinTest = ClinTest,
             YTrain = YTrain, YTest = YTest)

  saveRDS(res, paste0(dirData, 'splitted_', file))

  if(keep == TRUE){
    return(res)
  }
}


prepareData = function(FilePath, target = c('ages', 'AGE_YEARS', 'SUBSET_ANTIBIOTIC_HISTORY'), fs = FALSE){
  
  data = readRDS(FilePath)
  
  # Creating train and test data
  if (fs == TRUE) {
    featsFS = readRDS('~/git/Metagenomics/FS_lasso_genera.rds')
    featsFS = gsub('X', '', featsFS)
    train = cbind.data.frame(data$ClinTrain, data$OtuTrain[, intersect(featsFS, names(data$OtuTrain))], target = data$YTrain[, target])
  } else if(fs == FALSE){
    train = cbind.data.frame(data$ClinTrain, data$OtuTrain, target = data$YTrain[, target])
  }
  
  test = cbind.data.frame(data$ClinTest, data$OtuTest, target = data$YTest[, target])
  
  # Make names
  names(train) = make.names(names(train))
  names(test) = make.names(names(test))
  
  res = list(train = train, test = test)
}


preproRegr = function(data, target, ntile){

  require(dplyr)
  require(BurStMisc)

  # Convert ages into numeric
  if (ntile == TRUE){
    data$AGE_YEARS = ntile(data$AGE_YEARS, 60, result = 'numeric')
  }

  if (target == 'ages'){
    data$ages = as.numeric(data$ages)
  }

  # Change cvrts format  
  if ('ALCOHOL_CONSUMPTION' %in% names(data)) {
    data = mutate(data, ALCOHOL_CONSUMPTION = ifelse(data$ALCOHOL_CONSUMPTION == TRUE, 'yes', 'no'))
    data$ALCOHOL_CONSUMPTION = as.factor(data$ALCOHOL_CONSUMPTION)
  }
  
  if ('SUBSET_ANTIBIOTIC_HISTORY' %in% names(data)) {
    data = mutate(data, SUBSET_ANTIBIOTIC_HISTORY = ifelse(data$SUBSET_ANTIBIOTIC_HISTORY == TRUE, 'yes', 'no'))
    data$SUBSET_ANTIBIOTIC_HISTORY = as.factor(data$SUBSET_ANTIBIOTIC_HISTORY)
  }
  
  if ('SUBSET_HEALTHY' %in% names(data)) {
    data = mutate(data, SUBSET_HEALTHY = ifelse(data$SUBSET_HEALTHY == TRUE, 'yes', 'no'))
    data$SUBSET_HEALTHY = as.factor(data$SUBSET_HEALTHY)
  }
  
  if ('APPENDIX_REMOVED' %in% names(data)) {
    data = mutate(data, APPENDIX_REMOVED = ifelse(data$APPENDIX_REMOVED == TRUE, 'yes', 'no'))
    data$APPENDIX_REMOVED = as.factor(data$APPENDIX_REMOVED)
  }

  return(data)
}

preproClassif = function(data, target){

  require(dplyr)

  if (target == 'SUBSET_ANTIBIOTIC_HISTORY') {
    #data = subset(data, select = -c(SUBSET_ANTIBIOTIC_HISTORY, ANTIBIOTIC_HISTORY))
    data = data[, -c(8,9)]
  }
  
  if ('ALCOHOL_CONSUMPTION' %in% names(data)) {
    data = mutate(data, ALCOHOL_CONSUMPTION = ifelse(data$ALCOHOL_CONSUMPTION == TRUE, 'yes', 'no'))
    data$ALCOHOL_CONSUMPTION = as.factor(data$ALCOHOL_CONSUMPTION)
  }
  
  if ('SUBSET_ANTIBIOTIC_HISTORY' %in% names(data)) {
    data = mutate(data, SUBSET_ANTIBIOTIC_HISTORY = ifelse(data$SUBSET_ANTIBIOTIC_HISTORY == TRUE, 'yes', 'no'))
    data$SUBSET_ANTIBIOTIC_HISTORY = as.factor(data$SUBSET_ANTIBIOTIC_HISTORY)
  }
  
  if ('SUBSET_HEALTHY' %in% names(data)) {
    data = mutate(data, SUBSET_HEALTHY = ifelse(data$SUBSET_HEALTHY == TRUE, 'yes', 'no'))
    data$SUBSET_HEALTHY = as.factor(data$SUBSET_HEALTHY)
  }
  
  if ('APPENDIX_REMOVED' %in% names(data)) {
    data = mutate(data, APPENDIX_REMOVED = ifelse(data$APPENDIX_REMOVED == TRUE, 'yes', 'no'))
    data$APPENDIX_REMOVED = as.factor(data$APPENDIX_REMOVED)
  }
  
  # Delete order factors and change them for only factors
  if ('EXERCISE_FREQUENCY' %in% names(data)) {
    data$EXERCISE_FREQUENCY = factor(data$EXERCISE_FREQUENCY, ordered = F)
  }
  
  if ('ALCOHOL_FREQUENCY' %in% names(data)) {
    data$ALCOHOL_FREQUENCY = factor(data$ALCOHOL_FREQUENCY, ordered = F)
  }

  if ('SMOKING_FREQUENCY' %in% names(data)) {
    data$SMOKING_FREQUENCY = factor(data$SMOKING_FREQUENCY, ordered = F)
  }
  
  if ('SLEEP_DURATION' %in% names(data)) {
    data$SLEEP_DURATION = factor(data$SLEEP_DURATION, ordered = F)
  }

  return(data)

}



phyloseq_to_edgeR = function(physeq, group, method="RLE", ...){
  require("edgeR")
  require("phyloseq")
  # Enforce orientation.
  if( !taxa_are_rows(physeq) ){ physeq <- t(physeq) }
  x = as(otu_table(physeq), "matrix")
  # Add one to protect against overflow, log(0) issues.
  x = x + 1
  # Check `group` argument
  if( identical(all.equal(length(group), 1), TRUE) & nsamples(physeq) > 1 ){
    # Assume that group was a sample variable name (must be categorical)
    group = get_variable(physeq, group)
  }
  # Define gene annotations (`genes`) as tax_table
  taxonomy = tax_table(physeq, errorIfNULL=FALSE)
  if( !is.null(taxonomy) ){
    taxonomy = data.frame(as(taxonomy, "matrix"))
  } 
  
  # Now turn into a DGEList
  y = DGEList(counts=x, group=group, genes=taxonomy, remove.zeros = TRUE, ...)
  # Calculate the normalization factors
  z = calcNormFactors(y, method=method)
  # Check for division by zero inside `calcNormFactors`
  if( !all(is.finite(z$samples$norm.factors)) ){
    stop("Something wrong with edgeR::calcNormFactors on this data,
         non-finite $norm.factors, consider changing `method` argument")
  }
  # Estimate dispersions
  return(estimateTagwiseDisp(estimateCommonDisp(z)))
  }