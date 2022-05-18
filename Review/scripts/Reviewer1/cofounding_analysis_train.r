require(phyloseq)

# Load train data
# ====
train = readRDS('~/git/IBDpred/01_Preprocessing/data/train__IBD.rds')

# Extract clinical data
# ===
clin = as.data.frame(train@sam_data)

# Load train probabilities
# ===
require(mlr)
m = readRDS('~/git/IBDpred/02_Training/results/kruskal.test_40train__IBD_40_glmnet.rds')
res = as.data.frame(m)
model = getBMRModels(m)
best = model[[1]][[1]][[which.max(res$auc)]]

features = best$features

train.df = t(train@otu_table)
colnames(train.df) = paste0('X', colnames(train.df))
train.df = train.df[, features]
train.df = apply(train.df, 2, function(x) log2(x + 1))
train.df = cbind.data.frame(train.df, target = get_variable(train, 'IBD'))
train.df = makeClassifTask(data = train.df, target = 'target')
train.task = normalizeFeatures(train.df)

pred = predict(best, train.task, type = 'prob')
pROC::roc(pred$data$truth, pred$data$prob.NO, ci = T, plot = F)

# Add predictions to clinical data
clin$pred = pred$data$prob.NO

# Plot
require(ggpubr)
require(viridis)
toCompare = c('BMI', 'ALCOHOL_CONSUMPTION', 'LIVER_DISEASE', 'CANCER', 'DIET_TYPE',
              'DIABETES', 'COUNTRY', 'SEX', 'SMOKING_FREQUENCY', 'PROBIOTIC_FREQUENCY',
              'ANTIBIOTIC_HISTORY', 'AUTOIMMUNE', 'pred')
summary(clin[, toCompare])

# Ages
p1 = ggscatter(data = clin, x = 'AGE_YEARS', y = 'pred', add = 'reg.line', conf.int = T,
          add.params = list(color = viridis(3)[1], fill = 'lightgray')) +
  stat_cor(method = 'spearman')

# BMI
clin$BMI[which(clin$BMI > 40)] = NA
p2 = ggscatter(data = clin, x = 'BMI', y = 'pred', add = 'reg.line', conf.int = T,
          add.params = list(color = viridis(3)[1], fill = 'lightgray')) +
  stat_cor(method = 'spearman')

# Alcohol consumption
d = data.frame(
  alcohol = clin$ALCOHOL_CONSUMPTION,
  pred = clin$pred
)
d = d[which(!is.na(d$alcohol)),]
p3 = ggboxplot(d, x = 'alcohol', y = 'pred', fill = viridis(2)) + stat_compare_means(label.y.npc = 0.7)

# Antibiotic
d = data.frame(
  antibiotic = clin$SUBSET_ANTIBIOTIC_HISTORY,
  pred = clin$pred
)
d = d[which(!is.na(d$antibiotic)),]
p4 = ggboxplot(d, x = 'antibiotic', y = 'pred', fill = viridis(2)) + stat_compare_means(label.y.npc = 0.7)

# PROBIOTIC
d = data.frame(
  probiotic = as.character(clin$PROBIOTIC_FREQUENCY),
  pred = clin$pred
)
d = d[-which(d$probiotic == 'Unknown' | d$probiotic == 'Unspecified'),]
d$probiotic = factor(d$probiotic, ordered = T,
                     levels = c('Never', 'Rarely (a few times/month)', 'Occasionally (1-2 times/week)', 'Regularly (3-5 times/week)', 'Daily'),
                     labels = c('Never' = 'Never',
                                'Rarely (a few times/month)' = 'Rarely',
                                'Occasionally (1-2 times/week)' = 'Occasionally',
                                'Regularly (3-5 times/week)' = 'Regularly',
                                'Daily' = 'Daily'))

pp = ggboxplot(d, x = 'probiotic', y = 'pred', fill = viridis(5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(method = 'anova')

# Appendix
d = data.frame(
  appendix = clin$APPENDIX_REMOVED,
  pred = clin$pred
)
d = d[which(!is.na(d$appendix)),]
p5 = ggboxplot(d, x = 'appendix', y = 'pred', fill = viridis(2)) + stat_compare_means(label.y.npc = 0.7)

# Sex
d = data.frame(
  sex = clin$SEX,
  pred = clin$pred
)
d = d[-which(d$sex == 'Unknown' | d$sex == 'unspecified'),]
p6 = ggboxplot(d, x = 'sex', y = 'pred', fill = viridis(2)) + stat_compare_means(label.y.npc = 0.7)

ggpubr::ggarrange(p1, p2, p3, p4, p5, p6, pp, ncol = 4, nrow = 2)