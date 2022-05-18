# Figure 2a
# =====
require(forestplot)
require(viridis)

tabletext = cbind(
  c('Cohort', 'Gevers', 'Morgan'),
  c("Problem", 'CD', 'UC'),
  c('#Samples', '300', '66'),
  c('#Features', '12', '9'))

AUC.glmnet.fcbf = structure(list(
  mean  = c(NA, 0.444, 0.6142),
  lower = c(NA, 0.383, 0.4677),
  upper = c(NA,  0.505, 0.7607)),
  .Names = c("mean", "lower", "upper"),
  row.names = c(NA, -3L),
  class = "data.frame")

AUC.rf.fcbf = structure(list(
  mean  = c(NA, 0.477, 0.5778),
  lower = c(NA, 0.452, 0.4308),
  upper = c(NA,  0.504, 0.7249)),
  .Names = c("mean", "lower", "upper"),
  row.names = c(NA, -3L),
  class = "data.frame")

AUC.glmnet.k40 = structure(list(
  mean  = c(NA, 0.5478, 0.5532),
  lower = c(NA, 0.4824, 0.4069),
  upper = c(NA, 0.6133, 0.6994)),
  .Names = c("mean", "lower", "upper"),
  row.names = c(NA, -3L),
  class = "data.frame")

AUC.rf.k40 = structure(list(
  mean  = c(NA, 0.7588, 0.7391),
  lower = c(NA, 0.7042, 0.6082),
  upper = c(NA, 0.8134, 0.8699)),
  .Names = c("mean", "lower", "upper"),
  row.names = c(NA, -3L),
  class = "data.frame")

forestplot::forestplot(tabletext,
                       new_page = T,
                       hrzl_lines = gpar(col = "#444444"),
                       mean = cbind(AUC.glmnet.k40[, 'mean'], AUC.rf.k40[, 'mean']),
                       lower = cbind(AUC.glmnet.k40[, 'lower'], AUC.rf.k40[, 'lower']),
                       upper = cbind(AUC.glmnet.k40[, 'upper'], AUC.rf.k40[, 'upper']),
                       col=fpColors(box=viridis(2)),
                       line.margin = 0.5,
                       boxsize = .25,
                       clip = c(0.3, 1),
                       is.summary = c(T, F, F),
                       zero = 0.5,
                       xlab = 'AUC',
                       lty.ci = c(1,2),
                       xticks = c(.3, .4, .5, .6, .7, .8, .9),
                       grid = T,
                       legend = c("Glmnet", "Random Forest"))



# Figure 2b
# ====
morgan = readRDS('~/git/IBDpred/02_Training/results/retrain/validation/k40_morgan_uc.rds')$pred_glmnet
gevers = readRDS('~/git/IBDpred/02_Training/results/retrain/validation/k40_gevers_cd.rds')$pred_glmnet

require(mlr)
require(ROCR)
require(pROC)
ROCpred1 = asROCRPrediction(morgan)
ROCRperf1 = ROCR::performance(ROCpred1, "tpr", "fpr")

ROCpred2 = asROCRPrediction(gevers)
ROCRperf2 = ROCR::performance(ROCpred2, "tpr", "fpr")

roc1 = roc(response = ROCpred1@labels[[1]], predictor = ROCpred1@predictions[[1]], ci = T)
roc2 = roc(response = ROCpred2@labels[[1]], predictor = ROCpred2@predictions[[1]], ci = T)
roc.list1 = list(`UC (Morgan et al.)` = roc1, `CD (Gevers et al.)` = roc2)
ggroc(roc.list1) +
  theme_bw() +
  geom_abline(intercept = 1) +
  ggtitle('Glmnet model') +
  theme(legend.title = element_blank())


# 
# plot(ROCRperf1, col = 'blue')
# plot(ROCRperf2, col = 'red', add = T)
# legend('bottomright',
#        legend = c('Ulcerative Colitis', 'Crohn´s Disease'),
#        col = c('blue', 'red'))



# Figure 2c
# ===
morgan = readRDS('~/git/IBDpred/02_Training/results/retrain/validation/k40_morgan_uc.rds')$pred_rf
gevers = readRDS('~/git/IBDpred/02_Training/results/retrain/validation/k40_gevers_cd.rds')$pred_rf

require(mlr)
require(ROCR)
ROCpred1 = asROCRPrediction(morgan)
ROCRperf1 = ROCR::performance(ROCpred1, "tpr", "fpr")

ROCpred2 = asROCRPrediction(gevers)
ROCRperf2 = ROCR::performance(ROCpred2, "tpr", "fpr")


roc1 = roc(response = ROCpred1@labels[[1]], predictor = ROCpred1@predictions[[1]], ci = T)
roc2 = roc(response = ROCpred2@labels[[1]], predictor = ROCpred2@predictions[[1]], ci = T)
roc.list1 = list(`UC (Morgan et al.)` = roc1, `CD (Gevers et al.)` = roc2)
ggroc(roc.list1) +
  theme_bw() +
  geom_abline(intercept = 1) +
  ggtitle('Random Forest model') +
  theme(legend.title = element_blank())


# plot(ROCRperf1, col = 'blue')
# plot(ROCRperf2, col = 'red', add = T)
# legend('bottomright',
#        legend = c('Ulcerative Colitis', 'CrohnÂ´s Disease'),
#        col = c('blue', 'red'))



# Figure 2d
# ====
require(ggplot2)
require(reshape2)

dd = readxl::read_excel('~/git/IBDpred/Figures/utils/genus_availability.xlsx', col_names = T)

ggplot(melt(dd), aes(Genera, variable, fill = Phylo, alpha = value)) +
  geom_tile(colour = "gray50") +
  scale_alpha_identity(guide = "none") +
  coord_equal(expand = 0) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
