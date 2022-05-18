## NOTE!! Change working directory if necessary!
setwd('~/git/IBDpred')

# First step
# 1. Glomerate by Genus
# 2. Outliers analysis
# 3. Labeling and balancing data
require(phyloseq)
source('R/isolationForest.r')
source('R/processing_pipeline_functions.r')

# Load phyloseq object
filePath = '01_Preprocessing/data/phyloseq.rds'
taxLevel = 'Rank6'

pseq = readRDS(filePath)

print(paste0('Glomerating OTUs by ', taxLevel, ' ...'))
pseqGlom = tax_glom(pseq, taxLevel)
pats = isoForest(pseqGlom)

saveRDS(pseqGlom, file = '01_Preprocessing/data/pseqGlomRank6.rds')
saveRDS(pats, file = '01_Preprocessing/data/keepPats_after_isoForest.rds')

pseq2 = prune_samples(pats, pseqGlom)
pseq2 = subset_taxa(pseq2, (Rank6 != "g__"))

target = 'IBD'
res = labeling(pseq2, target)
saveRDS(res, file = paste0('01_Preprocessing/data/phyloseq_', target, '.rds'))


# =====
#  2º Step -- Split into train and test
dirData = '01_Preprocessing/data/'
setwd(dirData)

files = list.files(pattern = 'phyloseq_')
targets = gsub('.rds', '', substr(files, 10, 50))

for (i in seq_along(files)) {
  
  pseq = readRDS(files[i])
  smp_size = floor(0.85 * nsamples(pseq))
  set.seed(1993)
  train_ind = sample(seq_len(nsamples(pseq)), size = smp_size)
  trainPats = sample_names(pseq)[train_ind]
  testPats = sample_names(pseq)[-train_ind]
  
  train = subset_samples(pseq, sample_names(pseq) %in% trainPats)
  test = subset_samples(pseq, sample_names(pseq) %in% testPats)
  
  saveRDS(train, file = paste0(dirData, 'train_', targets[i], '.rds'))
  saveRDS(test, file = paste0(dirData, 'test_', targets[i], '.rds'))
  
  print(paste0(targets[i], ' done!'))
  
}


# =====
#  2º Step -- DEG analysis
dirData = '01_Preprocessing/data/'
setwd(dirData)

files = list.files(pattern = 'train_')
targets = gsub('.rds', '', substr(files, 7, 50))

for (i in seq_along(files)) {

  train = readRDS(files[i])
  degs = deg(train, targets[i])
  saveRDS(degs, file = paste0(dirData, 'deg_', targets[i], '.rds'))
  print(paste0(targets[i], ' done!'))
  
}


# =====
#  3º Step -- Add clinical variable and save them to run ML
dirData = '01_Preprocessing/data/'
setwd(dirData)
files = list.files(pattern = 'deg_')

for (i in seq_along(files)) {
  
  pseq = readRDS(files[i])
  otu = as.data.frame(t(otu_table(pseq)))
  
  otu = apply(otu, 2, function(x) log2(x + 1))
  
  target = gsub('.rds', '', substr(files[i], 5, 50))
  data = cbind.data.frame(otu, target = get_variable(pseq, target))

  toRun = '01_Preprocessing/data/toRun/'
  if (dir.exists(toRun) == FALSE) {
    dir.create(toRun)
    message('Creating toRun directory ...')
  }
  
  saveRDS(data, file = paste0(dirData, 'toRun/otu_', target))
  
}



# =====
#  4º Step -- Feature Selection
source('R/fs_ml_functions.r')
require(FCBF)

dirData = '01_Preprocessing/data/'
setwd(dirData)
outDir = '01_Preprocessing/data/toRun/'

files = list.files(pattern = 'train_')
targets = gsub('.rds', '', substr(files, 7, 50))

for (i in seq_along(files)) {
  
  data = readRDS(files[i])
  otu = as.data.frame(t(otu_table(data)))
  otu = apply(otu, 2, function(x) log2(x + 1))

  otu = cbind.data.frame(otu, target = get_variable(data, targets[i]))
  
  names(otu) = make.names(names(otu))
  
  fcbf = fast.cor.FS(otu, thresh = 0.0025)
  fs.abs(otu, 'kruskal.test', c(5, 10, 20, 40), outDir, gsub('.rds', '', files[i]))
  
  
  saveRDS(fcbf, file = paste0(outDir, 'fcbf_', targets[i]))
  print(paste0(targets[i], ' done!'))
  
}




