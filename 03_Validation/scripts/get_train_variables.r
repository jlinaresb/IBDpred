## NOTE!! Change working directory if necessary!
setwd('~/git/IBDpred')

require(mlr)

files = list.files('02_Training/results/', pattern = 'rf.rds')
trainVariables = list()
for (i in seq_along(files)) {
  r = readRDS(paste0('02_Training/results/', files[i]))
  m = getBMRModels(r)
  trainVariables[[i]] = m[[1]][[1]][[1]]$features
}


names(trainVariables) = sapply(strsplit(files, '_'), '[[', 1)
names(trainVariables) = c('fcbf', 'ldm', 'k10', 'k20', 'k40', 'k5', 'deg')
saveRDS(trainVariables, file = '03_Validation/data/trainVariables.rds')

