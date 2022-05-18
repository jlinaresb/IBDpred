## NOTE!! Change working directory if necessary!
setwd('~/git/IBDpred')

# Arguments
input.dir.path = '01_Preprocessing/data/toRun/'
path.algs = '02_Training/scripts/MachineLearning/'
pattern.algs = '_train.r'

ExperimentName = '11012021_genera_v3'

part = 'shared'
qos = 'shared_short'
time = '02:00:00'
nodes = 1
ntasks = 20

out.path = paste('02_Training/results/', ExperimentName, '/', sep = '')
out.filename.glmnet='glmnet.rds'
out.filename.rf='rf.rds'

if (dir.exists(out.path) == FALSE) {
  dir.create(out.path)
  message(paste('Creating ', ExperimentName, ' directory!'))
}


# Glmnet
lambda = c(0.0001,0.001,0.01,0.1,1)
alpha = c(0,0.15,0.25,0.35,0.5,0.65,0.75,0.85,1)
# Random Forest
mtry = c(2:8)
ntree = 1000
nodesize = c(1:3)

cv.inner = 'Holdout'
cv.outer = 'RepCV'

predict = c('both') # train or both
iters = 10
reps = 5
folds = 5
strat= TRUE
