# Prediction in Morgan
# ====

require(sva)
require(mlr)
require(ROCR)

# Arguments
# ==================
cohort = 'gevers' #gevers #morgan
ibd = 'cd' #uc #cd
fs = 'k40' #fcbf #k5 #k10 #k20 #k40 #deg #ldm

# Retrieve BMR models and validation data
# ===================
setwd('~/git/IBDpred/02_Training/results/retrain/')
files = list.files()
files = files[grep(fs, files)]
files = files[grep(ibd, files)]
files = files[grep(cohort, files)]
gl = readRDS(files[grep('glmnet', files)])
rf = readRDS(files[grep('rf', files)])

m = getBMRModels(rf)
bb = which.max(as.data.frame(rf)$acc)
best.rf = m[[1]][[1]][[bb]]

m = getBMRModels(gl)
bb = which.max(as.data.frame(gl)$acc)
best.gl = m[[1]][[1]][[bb]]

val = readRDS(paste0('~/projects/Metagenomics/validation/validation_IBD_', cohort, '_', ibd, '_', fs)) ; print(length(names(val)) - 1)
val = makeClassifTask(data = val, target = 'target')
val = normalizeFeatures(val)

# performance(predict(best.gl, val), measures = list(auc))
pred.gl = predict(best.gl, val, type = 'prob')
pROC::roc(pred.gl$data$truth, pred.gl$data$prob.NO, ci = T, plot = F)

# performance(predict(best.rf, val), measures = list(auc))
pred.rf = predict(best.rf, val, type = 'prob')
pROC::roc(pred.rf$data$truth, pred.rf$data$prob.NO, ci = T, plot = F)


saveRDS(list(pred_glmnet = pred.gl,
             pred_rf = pred.rf),
        file = paste0('~/git/IBDpred/02_Training/results/retrain/validation/', fs, '_', cohort, '_', ibd, '.rds'))