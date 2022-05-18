# See validation in train (merge = morgan + gevers)
# ===
# Not run!!
# setwd('~/git/IBDpred/Review/results/')
# 
# require(mlr)
# 
# files = list.files()
# bmrs = list()
# for (i in seq_along(files)) {
#   bmrs[[i]] = readRDS(files[i])
# }
# bmr = mergeBenchmarkResults(bmrs)
# saveRDS(bmr, file = 'merge_bmr_results.rds')


# Predict in AGP
# ====
# validation data
agp = readRDS('~/git/IBDpred/Review/data/agp_df.rds')

# bmr models
bmr = readRDS('~/git/IBDpred/Review/results/merge_bmr_results.rds')

plotBMRSummary(bmr, measure = auc, pretty.names = T)
plotBMRBoxplots(bmr, style = 'violin', measure = auc)

# select best model
models = getBMRModels(bmr)
folds = as.data.frame(bmr)
rf.folds = folds[which(folds$learner.id == 'classif.randomForest.tuned'),]
glmnet.folds = folds[which(folds$learner.id == 'classif.glmnet.tuned'),]

i.rf = as.data.frame(rf.folds)[which.max(as.data.frame(rf.folds)$auc),]
i.glmnet = as.data.frame(glmnet.folds)[which.max(as.data.frame(glmnet.folds)$auc),]

best.rf = models$kruskal.test_20kruskal_merge_nfeat_20$classif.randomForest.tuned[[i.rf$iter]]
best.glmnet = models$kruskal.test_40kruskal_merge_nfeat_40$classif.glmnet.tuned[[i.glmnet$iter]]

val = agp[, match(c(best.rf$features, 'target'), colnames(agp))]
val = makeClassifTask(data = val, target = 'target')
val = normalizeFeatures(val)
pred = predict(best.rf, val, type = 'prob')
pROC::roc(pred$data$truth, pred$data$prob.NO, ci = T, plot = T)

val = agp[, match(c(best.glmnet$features, 'target'), colnames(agp))]
val = makeClassifTask(data = val, target = 'target')
val = normalizeFeatures(val)
pred = predict(best.glmnet, val, type = 'prob')
pROC::roc(pred$data$truth, pred$data$prob.NO, ci = T, plot = T)
