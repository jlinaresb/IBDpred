# Functions to FS and ML

# FEATURE SELECTION Functions
# ===========================

# FCBF
fast.cor.FS = function(data, thresh){
  
  stopifnot('target' %in% names(data))
  # require(FCBF)
  
  y = as.factor(data$target)
  x = subset(data, select = - c(target))
  
  dis = discretize_exprs(t(x))
  # su_plot(dis, y)
  
  fcbf = fcbf(dis, y, verbose = T, thresh)
  xx = x[,fcbf$index]
  
  xx = as.data.frame(cbind(xx, target = y))
  
  return(xx)
  
}

fs.abs = function(data, fs.type, nfeat, output.dir, filename){
  
  stopifnot('target' %in% names(data))
  load('~/git/r-package-tcga-methodology/3.Filter/Mathematical/filter.methods.RData')
  stopifnot(fs.type %in% filter.methods$id)
  
  require(mlr)
  task = makeClassifTask(data = data, target = 'target')
  
  tasks = lapply(nfeat, function(x) filterFeatures(task, method = fs.type, abs = x))
  
  for (i in 1:length(nfeat)) {
    tasks[[i]]$task.desc$id =  paste(fs.type, ncol(tasks[[i]]$env$data) - 1 , sep = "_")
  }
  
  t = list()
  for (i in 1:length(tasks)) {
    t[[i]] = tasks[[i]]$env$data
    names(t)[[i]] = paste(fs.type, nfeat[i], sep = '_')
  }
  
  for (i in 1:length(t)) {
    saveRDS(t[[i]], file = paste(output.dir, names(t)[i], filename, sep = ''))
  }
  
  return(t)
}

# Spearman correlation
spearman = function(otu, cvrts, target, nfeat, outDir, expName){

	s = apply(otu, 2, function(x) cor.test(x, target, method = 'spearman')$estimate)
	s = sort(abs(s), decreasing = T)
	n = names(s)

	newDir = paste0(outDir, expName, '/')
	if (dir.exists(newDir) == FALSE){
		dir.create(newDir)
	}

	for (i in seq_along(nfeat)) {
		feats = n[seq(nfeat[i])]
		res = otu[, feats]
		res = cbind.data.frame(res, target = y)

		if (is.null(cvrts) == FALSE){
			res = cbind.data.frame(cvrts, res)
		}

		saveRDS(res, file = paste0(newDir, expName, '_', nfeat[i], '.rds'))
	}

	return(res)
}


filtering = function(data, fs.type, nfeat, Yvar, outDir, cvrts = FALSE, covariates){

  require(mlr)
  stopifnot('target' %in% names(data))

  task = makeClassifTask(data = data, target = 'target')

  tasks = lapply(nfeat, function(x) filterFeatures(task, method = fs.type, abs = x))

  for (i in 1:length(nfeat)) {
    tasks[[i]]$task.desc$id =  paste(fs.type, ncol(tasks[[i]]$env$data) - 1 , sep = "_")
  }

  t = list()
  for (i in 1:length(tasks)) {

    if (cvrts == TRUE){
      t[[i]] = tasks[[i]]$env$data
      names(t)[[i]] = paste(fs.type, nfeat[i], 'cvrts', sep = '_')
      res = cbind.data.frame(covariates, t[[i]])
      saveRDS(res,file = paste0(outDir, Yvar, '_', names(t)[i]))
    } else{
      t[[i]] = tasks[[i]]$env$data
      names(t)[[i]] = paste(fs.type, nfeat[i], sep = '_')
      saveRDS(t[[i]], file = paste0(outDir, Yvar, '_', names(t)[i]))
    }
  }
}





# MACHINE LEARNING Functions
# ===========================

# glmnet regression
RegrGlmnet = function(path, Yvar = c('AGE_YEARS'), outDir, keep = FALSE, save = TRUE, parallel=TRUE){

  require(mlr)
  library(parallelMap)
  
  name = basename(gsub('.rds', '', path))
  data = readRDS(path)  
  
  # Convert characters into factors
  ch = sapply(data, class) == 'character'
  it = which(ch == TRUE)

  if (length(it) > 0){
    for (i in seq_along(it)) {
      data[, it[i]] = as.factor(data[, it[i]])
    }
  }

  names(data) = make.names(names(data))
  
  # Create Task
  task = makeRegrTask(data = data, target = 'target')
  print('Removing Constant Features')
  task = removeConstantFeatures(task)
  print('Normalizing Features')
  task = normalizeFeatures(task)
  print('Select Hyperparameters')
  
  # Hyperparameter tuning
  ctrl<-makeTuneControlGrid()
  inner = makeResampleDesc(method = 'Holdout', predict = 'both')
  outer = makeResampleDesc(method = 'RepCV', predict = 'both', reps = 5, folds = 10)
  
  lambda = c(0.0001,0.001,0.01,0.1,1)
  alpha = c(0,0.15,0.25,0.35,0.5,0.65,0.75,0.85,1)
  
  psglmnet = makeParamSet(
    makeDiscreteParam("lambda", lambda),
    makeDiscreteParam("alpha",alpha)
  )
  l<-makeLearner("regr.glmnet")
  lrn_glmnet<-makeTuneWrapper(l, inner, psglmnet, measures = rrse, ctrl, show.info=T)
  
  print('Training the model')
  
  if(parallel == TRUE){
    parallelStartMulticore(20L , level = 'mlr.tuneParams')
  }
  
  bmr = benchmark(lrn_glmnet, task , outer , measures =  list(rmse, rrse, rsq, mse, kendalltau, expvar, mae, rae, spearmanrho) , show.info = T , models = T)
  parallelStop()

  if (keep == TRUE) {
    return(list(trainData = task, res = bmr))
  }
  
  if (save == TRUE) {
    if (dir.exists(outDir) == FALSE) {
      dir.create(outDir)
      message('Creating directory --Exec-- ...')
    }
    saveRDS(list(trainData = task, res = bmr), file = paste0(outDir, '/glmnet_', Yvar, name, '_fs', '.rds'))
  }
}


# Random Forest regression
RegrRF = function(path, Yvar = c('AGE_YEARS'), outDir, keep = FALSE, save = TRUE, parallel=TRUE){

  require(mlr)
  library(parallelMap)
  
  name = basename(gsub('.rds', '', path))
  data = readRDS(path) 

  # Convert characters into factors
  ch = sapply(data, class) == 'character'
  it = which(ch == TRUE)

  if (length(it) > 0){
    for (i in seq_along(it)) {
      data[, it[i]] = as.factor(data[, it[i]])
    }
  }

  names(data) = make.names(names(data))
  
  # Create Task
  task = makeRegrTask(data = data, target = 'target')
  print('Removing Constant Features')
  task = removeConstantFeatures(task)
  print('Normalizing Features')
  task = normalizeFeatures(task)
  print('Select Hyperparameters')
  
  # Hyperparameter tuning
  ctrl<-makeTuneControlGrid()
  inner = makeResampleDesc(method = 'Holdout', predict = 'both')
  outer = makeResampleDesc(method = 'RepCV', predict = 'both', reps = 5, folds = 10)

  if (ncol(data) < 15){
    if (ncol(data) <= 4){
      mtry = c(2,3)
    } else{
      mtry = c(2:4)
    }
  } else{
    mtry = round(seq(sqrt(ncol(data)) - 3, sqrt(ncol(data)) + 3))
  }
  
  ntree = 1000
  nodesize = c(1:5)
  
  psrf = makeParamSet(
      makeDiscreteParam("mtry", values = mtry),
      makeDiscreteParam("ntree", values= ntree),
      makeDiscreteParam("nodesize", values= nodesize))

  l = makeLearner("regr.randomForest")
  lrn_rf = makeTuneWrapper(l, inner, psrf, measures = rrse, ctrl, show.info=T)
  
  print('Training the model')

  if(parallel == TRUE){
    parallelStartMulticore(20L , level = 'mlr.tuneParams')
  }
  
  bmr = benchmark(lrn_rf, task , outer , measures =  list(rmse, rrse, rsq, mse, kendalltau, expvar, mae, rae, spearmanrho) , show.info = T , models = T)
  parallelStop()

  if (keep == TRUE) {
    return(list(trainData = task, res = bmr))
  }
  
  if (save == TRUE) {
    if (dir.exists(outDir) == FALSE) {
      dir.create(outDir)
      message('Creating directory --Exec-- ...')
    }
    saveRDS(list(trainData = task, res = bmr), file = paste0(outDir, '/rf_', Yvar, name, '.rds'))
  }
}

# SVM regression
RegrSVM = function(path, Yvar = c('AGE_YEARS'), outDir, keep = FALSE, save = TRUE, parallel=TRUE){

  require(mlr)
  library(parallelMap)

  name = basename(gsub('.rds', '', path))
  data = readRDS(path) 

  # Convert characters into factors
  ch = sapply(data, class) == 'character'
  it = which(ch == TRUE)

  if (length(it) > 0){
    for (i in seq_along(it)) {
      data[, it[i]] = as.factor(data[, it[i]])
    }
  }

  if ('ALCOHOL_FREQUENCY' %in% names(data)){
    data$ALCOHOL_FREQUENCY = factor(data$ALCOHOL_FREQUENCY, ordered = F)
  }

  names(data) = make.names(names(data))  
  
  # Create Task
  task = makeRegrTask(data = data, target = 'target')
  print('Removing Constant Features')
  task = removeConstantFeatures(task)
  print('Normalizing Features')
  task = normalizeFeatures(task)
  print('Select Hyperparameters')
  
  # Hyperparameter tuning
  ctrl<-makeTuneControlGrid()
  inner = makeResampleDesc(method = 'Holdout', predict = 'both')
  outer = makeResampleDesc(method = 'RepCV', predict = 'both', reps = 5, folds = 10)
  
  C = 2^c(-12:12)
  sigma = 2^c(-12:12)
  
  # Hyperparameter Tuning
  pssvm = makeParamSet(
    makeDiscreteParam("C", values = C),  
    makeDiscreteParam("sigma", values= sigma)
  )
  l = makeLearner("regr.ksvm")
  lrn_ksvm = makeTuneWrapper(l , inner, pssvm, measures = rrse, ctrl, show.info = T)

  print('Training the model')

  if(parallel == TRUE){
    parallelStartMulticore(20L , level = 'mlr.tuneParams')
  }
  
  bmr = benchmark(lrn_ksvm, task , outer , measures =  list(rmse, rrse, rsq, mse, kendalltau, expvar, mae, rae, spearmanrho) , show.info = T , models = T)
  parallelStop()

  if (keep == TRUE) {
    return(list(trainData = task, res = bmr))
  }
  
  if (save == TRUE) {
    if (dir.exists(outDir) == FALSE) {
      dir.create(outDir)
      message('Creating directory --Exec-- ...')
    }
    saveRDS(list(trainData = task, res = bmr), file = paste0(outDir, '/svm_', Yvar, name, '.rds'))
  }
}


# XGBoost regression
RegrXGB = function(path, Yvar = c('AGE_YEARS'), outDir, keep = FALSE, save = TRUE, parallel=TRUE){

  require(mlr)
  library(parallelMap)
  
  name = basename(gsub('.rds', '', path))
  data = readRDS(path)  

  # Convert variables according algorithms specifications
  ch = sapply(data, class) == 'character'
  it = which(ch == TRUE)

  if (length(it) > 0){
    for (i in seq_along(it)) {
      data[, it[i]] = as.factor(data[, it[i]])
      data[, it[i]] = as.vector(as.numeric(data[, it[i]]))
    }
  }

  if ('ALCOHOL_FREQUENCY' %in% names(data)){
    data$ALCOHOL_FREQUENCY = as.vector(as.numeric(data$ALCOHOL_FREQUENCY))
  }
  
  if ('DIET_TYPE' %in% names(data)){
    data$DIET_TYPE = as.vector(as.numeric(data$DIET_TYPE))
  }  

  names(data) = make.names(names(data))
  
  # Create Task
  task = makeRegrTask(data = data, target = 'target')
  print('Removing Constant Features')
  task = removeConstantFeatures(task)
  print('Normalizing Features')
  task = normalizeFeatures(task)
  print('Select Hyperparameters')
  
  # Hyperparameter tuning
  ctrl<-makeTuneControlGrid()
  inner = makeResampleDesc(method = 'Holdout', predict = 'both')
  outer = makeResampleDesc(method = 'RepCV', predict = 'both', reps = 5, folds = 10)
 
  # Hyperparameter Tuning
  psxgb = makeParamSet(
    # The number of trees in the model (each one built sequentially)
    makeIntegerParam("nrounds", lower = 100, upper = 500),
    # number of splits in each tree
    makeIntegerParam("max_depth", lower = 1, upper = 10),
    # "shrinkage" - prevents overfitting
    makeNumericParam("eta", lower = .1, upper = .5),
    # L2 regularization - prevents overfitting
    makeNumericParam("lambda", lower = -1, upper = 0, trafo = function(x) 10^x)
  )
  l = makeLearner("regr.xgboost")
  lrn_xgb = makeTuneWrapper(l , inner, psxgb, measures = rrse, ctrl, show.info = T)

  print('Training the model')

  if(parallel == TRUE){
    parallelStartMulticore(20L , level = 'mlr.tuneParams')
  }
  
  bmr = benchmark(lrn_xgb, task , outer , measures =  list(rmse, rrse, rsq, mse, kendalltau, expvar, mae, rae, spearmanrho) , show.info = T , models = T)
  parallelStop()

  if (keep == TRUE) {
    return(list(trainData = task, res = bmr))
  }
  
  if (save == TRUE) {
    if (dir.exists(outDir) == FALSE) {
      dir.create(outDir)
      message('Creating directory --Exec-- ...')
    }
    saveRDS(list(trainData = task, res = bmr), file = paste0(outDir, '/xgboost_', Yvar, name, '.rds'))
  }
}


# Glmnet classification
ClassifGlmnet = function(path, target = c('ages', 'SUBSET_ANTIBIOTIC_HISTORY'), outDir, keep = TRUE, save = FALSE){
  
  require(mlr)
  
  name = strsplit(basename(path), '_')[[1]][4]
  
  # Preparing data (load train and test data)
  data = prepareData(path, target)
  data = data$train
  names(data) = make.names(names(data))
  
  data = preproClassif(data, target)
  
  print('Making task')
  task =  makeClassifTask(data = data, target = target)
  
  print('Removing Constant Features')
  task = removeConstantFeatures(task)
  
  print('Normalizing Features')
  task = normalizeFeatures(task)
  
  print('Select Hyperparameters')
  print(task)
  
  # Estrategia de búsqueda
  ctrl<-makeTuneControlGrid()  
  inner = makeResampleDesc(method = 'Holdout', predict = 'both', stratify = T)
  
  lambda = c(0.0001,0.001,0.01,0.1,1)
  alpha = c(0,0.15,0.25,0.35,0.5,0.65,0.75,0.85,1)
  
  # Hyperparameter Tuning
  psglmnet = makeParamSet(
    makeDiscreteParam("lambda", lambda),
    makeDiscreteParam("alpha",alpha)
  )
  l<-makeLearner("classif.glmnet", predict.type = "prob")
  lrn_glmnet<-makeTuneWrapper(l, inner, psglmnet, measures = acc, ctrl, show.info=T)
  
  print('Select outer resampling')
  
  # Outer Cross-Validation
  outer = makeResampleDesc(method = 'RepCV', predict = 'both', reps = 5, folds = 10, stratify = T)
  
  print('Training the model')
  
  if (!('parallelMap' %in% installed.packages()[,"Package"])){
    message('Installing packages...')
    install.packages('parallelMap')
    library(parallelMap)
    parallelStartMulticore(20L , level = 'mlr.tuneParams')
  } else {
    library(parallelMap)
    parallelStartMulticore(20L , level = 'mlr.tuneParams')
    # parallelStartSocket(6, level = 'mlr.tuneParams')
  }
  
  # Benchmarking
  bmr = benchmark(lrn_glmnet, task , outer , measures =  list(auc, acc , mmce) , show.info = T , models = T)
  parallelStop()

  if (keep == TRUE) {
    return(list(trainData = task, res = bmr))
  }
  
  if (save == TRUE) {
    saveRDS(list(trainData = task, res = bmr), file = paste0(outDir, 'glmnetClassif_', target, name, '.rds'))
  }
}