rf.bmr.slurm = function(data, name, path = '', filename = '', cv.inner, cv.outer, mtry, ntree, nodesize){
  
  require(mlr)
  print('Making task')
  task =  makeClassifTask(id = paste(name, 'nfeat', ncol(data)-1, sep = '_'), data = data, target = 'target')
  
  print('Removing Constant Features')
  task = removeConstantFeatures(task)
  
  print('Normalizing Features')
  #task = normalizeFeatures(task)
  
  print('Select Hyperparameters')
  
  # Estrategia de búsqueda
  ctrl<-makeTuneControlGrid()
  
  if (cv.inner == 'CV') {
    inner = makeResampleDesc(method = 'CV', predict = 'both', iters = 10, stratify = TRUE)
  } else if (cv.inner == 'RepCV'){
    inner = makeResampleDesc(method = 'RepCV', predict = 'both', reps = 5, folds = 10, stratify = TRUE)
  } else{
    inner = makeResampleDesc(method = 'Holdout', predict = 'both', stratify = TRUE)
  }
  
  # Hyperparameter Tuning
  psrf<-makeParamSet(
    #makeDiscreteParam("mtry", values = c((round(sqrt(ncol(task$env$data))) - mtry) : (round(sqrt(ncol(task$env$data))) + mtry))),
    makeDiscreteParam("mtry", values = mtry),
    makeDiscreteParam("ntree", values= ntree),
    makeDiscreteParam("nodesize", values= nodesize)
  )
  l<-makeLearner("classif.randomForest", predict.type = "prob")
  lrn_rf<-makeTuneWrapper(l,  resampling = inner, par.set = psrf, measures = auc, control=ctrl,  show.info = T)
  
  print('Select outer resampling')
  
  # Outer Cross-Validation
  if (cv.outer == 'CV') {
    outer = makeResampleDesc(method = 'CV', predict = 'both', iters = 10, stratify = TRUE)
  } else if (cv.outer == 'RepCV'){
    outer = makeResampleDesc(method = 'RepCV', predict = 'both', reps = 5, folds = 10, stratify = TRUE)
  } else if (cv.outer == 'Holdout'){
    outer = makeResampleDesc(method = 'Holdout', predict = 'both', stratify = TRUE)
  } else{
    outer = makeResampleDesc(method = 'LOO', predict = 'both')
  }
  
  print('Training the model')
  
  if (!('parallelMap' %in% installed.packages()[,"Package"])){
    message('Installing packages...')
    install.packages('parallelMap')
    library(parallelMap)
    parallelStartMulticore(20L , level = 'mlr.tuneParams')
  } else {
    library(parallelMap)
    parallelStartMulticore(20L , level = 'mlr.tuneParams')
  }
  
  # Benchmarking
  bmr = benchmark(lrn_rf, task , outer , measures =  list(acc , auc, mmce) , show.info = T , models = T)
  
  saveRDS(bmr, file = paste(out.path, name, '_', ncol(data)-1, '_', out.filename.rf, sep = ''))
  
  parallelStop()
}