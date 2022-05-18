# Utils Plots
# ==============

require(plyr)

convertIDs = function(ids, pseq){
  pp = prune_taxa(ids, pseq)
  res = pp@tax_table@.Data[, 'Rank6']
  
  return(res)
}

RF.vi = function(bestModel, pseq, name){
  
  feats = getFeatureImportance(bestModel[[1]])$res$variable
  
  sumImp = rep(0, length(feats))
  for (i in seq_along(bestModel)) {
    iters = getFeatureImportance(bestModel[[i]])$res$importance
    sumImp =+ sumImp + iters
  }
  
  VarImp = data.frame(features = feats,
                      Importance = sumImp)
  rownames(VarImp) = VarImp$features

  VarImp$features = convertIDs(gsub('X', '', VarImp$features), pseq)
  tt = prune_taxa(gsub('X', '', feats), pseq)
  tt = tt@tax_table[,1:2]; rownames(tt) = paste0('X', rownames(tt))
  tt = tt[match(rownames(VarImp), rownames(tt)),]
  stopifnot(rownames(tt) == rownames(VarImp))
  VarImp = cbind.data.frame(VarImp, Phylum = as.vector(tt[,'Rank2']))
  
  p = ggplot(data = VarImp, aes(x = Importance, y = reorder(features, Importance), fill = Phylum)) + 
            geom_bar(stat = "identity", width = 0.5)
  p = p + scale_x_continuous(name = paste0("VI RF ", name))
  p = p + scale_y_discrete(name = "Features")
  p = p + theme_minimal()
  p = p + theme(axis.title.x = element_text(color = 'black', size = 15),
                axis.title.y = element_text(color = 'black', size = 15))
  p = p + theme(axis.text.x = element_text(angle = 30)) + theme(legend.position = 'bottom', legend.box = "horizontal")
  
  return(p)
}

glmnet.vi = function(bestModel, pseq, name){
  
  glmnet = bestModel
  sumBetas = rep(0, length(glmnet[[1]]$features))
  for (i in seq_along(bestModel)) {
    learner.models = getLearnerModel(glmnet[[i]])
    sumBetas = sumBetas + as.vector((learner.models$learner.model$beta))
  }
  
  names(sumBetas) = glmnet[[1]]$features
  
  df = data.frame(features = names(sumBetas),
                  Importance = sumBetas)
  
  df$features = convertIDs(gsub('X', '', df$features), pseq)
  tt = prune_taxa(gsub('X', '', rownames(df)), pseq)
  tt = tt@tax_table[,1:2]; rownames(tt) = paste0('X', rownames(tt))
  tt = tt[match(rownames(df), rownames(tt)),]
  stopifnot(rownames(tt) == rownames(df))
  df = cbind.data.frame(df, Phylum = as.vector(tt[,'Rank2']))

  # Plotting
  p = ggplot(data = df, aes(x = Importance, y = reorder(features, -Importance), fill = Phylum)) + geom_bar(stat = "identity", width = 0.5)
  p = p + scale_x_continuous(name = paste0("VI glmnet ", name))
  p = p + scale_y_discrete(name = "Features")
  p = p + theme_minimal()
  p = p + theme(axis.title.x = element_text(color = 'black', size = 15),
                axis.title.y = element_text(color = 'black', size = 15))
  p = p + theme(axis.text.x = element_text(angle = 30)) + theme(legend.position = 'bottom', legend.box = "horizontal")
  
  return(p)
  
}

selectBest = function(bmr, pseq){
  
  d = as.data.frame(bmr)
  aggs = ddply(d, .(task.id, learner.id), summarise, auc=mean(auc), acc = mean(acc), mmce = mean(mmce))
  
  maxAUC = which.max(aggs$auc)
  minMMCE = which.min(aggs$mmce)
  
  if ((maxAUC == minMMCE) == TRUE) {
    best = maxAUC
  } else {
    A = aggs$auc[maxAUC] + aggs$acc[maxAUC]
    B = aggs$auc[minMMCE] + aggs$acc[minMMCE]
    s = which.max(c(A,B))
    if (s == 1) {
      best = maxAUC
    } else{
      best = minMMCE
    }
  }
  
  task = aggs$task.id[best]
  learner = aggs$learner.id[best]
  print(learner)
  
  m = getBMRModels(bmr)
  m = m[[task]][[learner]]
  
  
  if (learner == 'classif.randomForest.tuned') {
    RF.vi(m, pseq, name = task)
  }else{
    glmnet.vi(m, pseq, name = task)
  }
  
}




barPlot_FromVenn = function(problem){
  
  # Arguments
  #   problem: which problem do you want to plot?
  # Results
  #   ggplot bar plot with counts of appaerences of each genera in each FS

  # Packages
  require(plyr)
  require(dplyr)
  require(phyloseq)
  library(tidyverse)
  require(ggVennDiagram)
  require(ggplot2)
  require(mlr)
  require(ggpubr)
  require(grid)
  require(ggupset)
  require(ggalt)
  
  source('~/git/Metagenomics/utils_toPlot.r')
  
  if (problem == 'ANTIBIOTIC_HISTORY') {
    p = 1
    testName = 'ANTIBIOTIC'
  } else if (problem == 'ASD'){
    p = 2
    testName = 'ASD'
  } else if (problem == 'COUNTRY'){
    p = 3
    testName = 'COUNTRY'
  } else if (problem == 'DIET_TYPE'){
    p = 4
    testName = 'DIET TYPE'
  } else if (problem == 'IBD'){
    p = 5
    testName = 'IBD'
  } else if (problem == 'IBS'){
    p = 6
    testName = 'IBS'
  } else if (problem == 'LUNG_DISEASE'){
    p = 7
    testName = 'LUNG'
  }
  
  inputDir = '~/projects/Metagenomics/data/phyloseq/targeted4/'
  pathTrain = paste0(inputDir, 'train__', problem, '.rds')
  pathTest = paste0(inputDir, 'test__', problem, '.rds')
  
  # Retrieve variables in train
  trainVars = readRDS('~/projects/Metagenomics/trainVariables.rds')
  ldm_trVars = readRDS('~/projects/Metagenomics/ldm_trainVariables.rds')
  trainVars$ldm = ldm_trVars
  
  # Retrieve results in validation
  colNames = c(rep('FCBF', 2), rep('K5', 2), rep('K10', 2), rep('K20', 2), rep('K40', 2), rep('LDM', 2), rep('DEG', 2))
  msr = c('.Acc', '.AUC')
  colNames = paste0(colNames, msr)
  colNames = c('FS', 'Algorithm', colNames)
  t = readxl::read_excel('~/projects/Metagenomics/results/test_results.xlsx', sheet = 3, skip = 5, col_names = colNames)
  pp = t$FS[-which(is.na(t$FS))]
  t$FS = c(rbind(pp, pp))
  
  # Read train and test sets
  train = readRDS(pathTrain)
  
  x = convertIDs(gsub('X', '', trainVars$fcbf[[p]]), train)
  y = convertIDs(gsub('X', '', trainVars$kruskal40[[p]]), train)
  z = convertIDs(trainVars$deg[[p]], train)
  w = convertIDs(gsub('X', '',trainVars$ldm[[p]]), train)
  
  v = list(FCBF = x,
           kruskal = y,
           DEG = z,
           LDM = w)
  
  venn = ggVennDiagram(v)
  
  
  # 2.c Barplot with most representative features (genera/phyla)
  n4 = strsplit(venn$plot_env$data$text[4], ';')[[1]]
  n3 = c(strsplit(venn$plot_env$data$text[3], ';')[[1]],
         strsplit(venn$plot_env$data$text[5], ';')[[1]],
         strsplit(venn$plot_env$data$text[7], ';')[[1]],
         strsplit(venn$plot_env$data$text[11], ';')[[1]])
  n2 = c(strsplit(venn$plot_env$data$text[2], ';')[[1]],
         strsplit(venn$plot_env$data$text[6], ';')[[1]],
         strsplit(venn$plot_env$data$text[8], ';')[[1]],
         strsplit(venn$plot_env$data$text[10], ';')[[1]],
         strsplit(venn$plot_env$data$text[12], ';')[[1]],
         strsplit(venn$plot_env$data$text[14], ';')[[1]])
  
  n4 = gsub(' ', '', n4)
  n4 = gsub('\n', '', n4)
  n3 = gsub(' ', '', n3)
  n3 = gsub('\n', '', n3)
  n2 = gsub(' ', '', n2)
  n2 = gsub('\n', '', n2)
  
  d = data.frame(genera = c(n4, n3, n2),
                 counts = c(rep(4, length(n4)), rep(3, length(n3)), rep(2, length(n2))))
  p2 = ggplot(data = d, aes(x = counts, y = reorder(genera, +counts))) +
    geom_bar(stat = 'identity', width = 0.3) + theme_minimal() +
    ylab(element_blank()) + xlab('Appearences in FS signatures')
  
  return(p2)
}

