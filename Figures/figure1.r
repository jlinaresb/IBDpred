# Arguments
# ===
setwd('~/git/IBDpred/')

problem = 'IBD'
cesga = F
saved = F
p = 5
testName = 'IBD'

# Packages
require(plyr)
require(dplyr)
require(phyloseq)
library(tidyverse)
require(ggplot2)
require(mlr)
require(ggpubr)
require(grid)
require(ggupset)
require(ggalt)

source('R/utils_to_plot.r')

inputDir = '01_Preprocessing/data/'
trainVarsPath = '03_Validation/data/trainVariables.rds'
tPath = '02_Training/results/test/test_results.xlsx'
bmrDirPath = '02_Training/results/mergedBMR/'
outDir = 'Figures/plots/'

pathTrain = paste0(inputDir, 'train__', problem, '.rds')
pathTest = paste0(inputDir, 'test__', problem, '.rds')

# Retrieve variables in train
trainVars = readRDS(trainVarsPath)

# Retrieve results in validation
colNames = c(rep('FCBF', 2), rep('K5', 2), rep('K10', 2), rep('K20', 2), rep('K40', 2), rep('LDM', 2), rep('DEG', 2))
msr = c('.Acc', '.AUC')
colNames = paste0(colNames, msr)
colNames = c('FS', 'Algorithm', colNames)
t = readxl::read_excel(tPath, sheet = 3, skip = 5, col_names = colNames)
pp = t$FS[-which(is.na(t$FS))]
t$FS = c(rbind(pp, pp))

# Read train and test sets
train = readRDS(pathTrain)
test = readRDS(pathTest)

# =============
# Plotting
# =============

# 2.a Table with data summary of train/test
t1 = data.frame(
  Set = c('Train', 'Test'),
  Taxa = c(ntaxa(train), ntaxa(test)),
  Samples = c(nsamples(train), nsamples(test)),
  Healthy = c(table(get_variable(train, problem))[1], table(get_variable(test, problem))[1]),
  IBD = c(table(get_variable(train, problem))[2], table(get_variable(test, problem))[2]))
t1 = ggtexttable(t1, rows = NULL)


# 2.b Venn Plot with features intersect
x = convertIDs(gsub('X', '', trainVars$fcbf), train)
y = convertIDs(gsub('X', '', trainVars$k40), train)
z = convertIDs(gsub('X', '', trainVars$deg), train)
w = convertIDs(gsub('X', '',trainVars$ldm), train)

v = list(FCBF = x,
         kruskal = y,
         DEG = z,
         LDM = w)

kk = transform(stack(v),name=substr(ind,1,1))
kk = group_by(.data = kk, values)

k1 =  summarize(.data = kk, FS = list(ind))
p1.1 = ggplot(data = k1, aes(x = FS)) + 
  geom_bar(position = 'dodge', width = 0.5) + 
  geom_text(position = position_dodge(width = 0.9),
    hjust = 0.5, vjust = 1, stat = 'count', colour = "#FFFFFF",
    aes(label=after_stat(count))) +
  scale_y_continuous(breaks = NULL,  name = "") + 
  scale_x_upset() +
  theme_minimal()

k2 = kk[, 2]
k2 = k2 %>% 
  select(ind) %>%
  unnest(cols = ind) %>%
  count(ind) %>% 
  mutate(ind = fct_reorder(as.factor(ind), n))

p1.2 = ggplot(data = k2, aes(x = ind, y = n, label = n)) +
  geom_col(position = 'dodge', width = 0.5) +
  geom_text(position = position_dodge(width = 0.9), hjust = 0, colour = "#FFFFFF") +
  coord_flip() +
  scale_y_reverse() +
  xlab("") + ylab("") +
  theme_void() +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
  
p1 = cowplot::plot_grid(
  cowplot::plot_grid(NULL, p1.2 + theme(plot.margin = unit(c(-2, 10, -5, 1), "pt")), ncol = 1, rel_heights = c(1.25, 0.5)),
  p1.1, nrow = 1, rel_widths = c(0.75,2)
)
 

# 2.d BMR Summary plot in train set
bmrDir = bmrDirPath
files = list.files(bmrDir)
bmr = readRDS(paste0(bmrDir, files[grep(problem, files)]))

d = as.data.frame(bmr)
aggs = ddply(d, .(task.id, learner.id), summarise, auc=mean(auc), acc = mean(acc), mmce = mean(mmce))
aggs$task.id[grep('fcbf', aggs$task.id)] = 'FCBF'
aggs$task.id[grep('fit', aggs$task.id)] = 'LDM'
aggs$task.id[grep('test_5', aggs$task.id)] = 'K5'
aggs$task.id[grep('test_10', aggs$task.id)] = 'K10'
aggs$task.id[grep('test_20', aggs$task.id)] = 'K20'
aggs$task.id[grep('test_40', aggs$task.id)] = 'K40'
aggs$task.id[grep('otu', aggs$task.id)] = 'DEG'

aggs$learner.id[grep('glmnet', aggs$learner.id)] = 'Glmnet'
aggs$learner.id[grep('randomForest', aggs$learner.id)] = 'RF'

aggs$task.id = factor(aggs$task.id, levels = rev(c('DEG', 'LDM', 'FCBF', 'K5', 'K10', 'K20', 'K40')))
names(aggs)[2] = 'Algorithm'

# p3 = ggplot(aggs, aes(x = auc, y = task.id, color = Algorithm)) +
#  geom_point(size = 7) + xlab('AUC') + ylab('Feature Selection') + xlim(0, 1)


# 2.e Table with data summary in test
test = t[which(t$FS == testName),]
tt = data.frame(FS = c(rep('FCBF', 2), rep('K5', 2), rep('K10', 2), rep('K20', 2), rep('K40', 2), rep('LDM', 2), rep('DEG', 2)),
                Algorithm = rep(c('RF', 'Glmnet'), 7),
                AUC = c(test$FCBF.AUC, test$K5.AUC, test$K10.AUC, test$K20.AUC, test$K40.AUC, test$LDM.AUC, test$DEG.AUC),
                Acc = c(test$FCBF.Acc, test$K5.Acc, test$K10.Acc, test$K20.Acc, test$K40.Acc, test$LDM.Acc, test$DEG.Acc))

tt$AUC = as.numeric(tt$AUC)
tt$Acc = as.numeric(tt$Acc)
tt$FS = factor(tt$FS, levels = rev(c('DEG', 'LDM', 'FCBF', 'K5', 'K10', 'K20', 'K40')))

#p4 = ggplot(tt, aes(x = AUC, y = FS, color = Algorithm)) + 
#  geom_point(size = 7) + xlim(0, 1) + ylab('Feature Selection')


dtrain = as.data.frame(aggs)
dtest = as.data.frame(tt)

dtrain = dtrain[order(dtrain$task.id),]
dtrain = dtrain[order(dtrain$Algorithm),]
dtest = dtest[order(dtest$FS),]
dtest = dtest[order(dtest$Algorithm),]

df = data.frame(fs.technique = dtrain$task.id,
                train = dtrain$auc,
                test = dtest$AUC,
                ml.model = dtrain$Algorithm)

p2 = ggplot(df, aes(y=fs.technique, x=train, xend=test, colour=ml.model)) +
  geom_dumbbell(size=3, colour_x = "#440154FF", colour_xend = "#FDE725FF", show.legend = T,
                dot_guide=TRUE, dot_guide_size=0.25, position=position_dodge2()) +
  labs(x=NULL, y=NULL, title=NULL,subtitle=NULL) +
  theme_minimal() +
  theme(panel.grid.major.x=element_line(size=0.05),legend.title = element_blank())



# 2.f Varible Importance of best model (50 iterations)
p3 = selectBest(bmr, train)

# ===========

if (saved == T){
  saveRDS(t1, file = paste0(outDir, 'panel_t1_', problem))
  saveRDS(p1, file = paste0(outDir, 'panel_p1_', problem))
  saveRDS(p2, file = paste0(outDir, 'panel_p2_', problem))
  saveRDS(p3, file = paste0(outDir, 'panel_p3_', problem))
} else {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow = 8, ncol = 2)))
  define_region <- function(row, col){
    viewport(layout.pos.row = row, layout.pos.col = col)
  } 
  
  print(t1, vp = define_region(row = 1, col = 1))
  print(p1, vp = define_region(row = 1:3, col = 2))
  print(p2, vp = define_region(row = 4:8, col = 2))
  print(p3, vp = define_region(row = 2:8, col = 1))
}
  


