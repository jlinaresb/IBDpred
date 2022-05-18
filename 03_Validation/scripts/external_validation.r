## NOTE!! Change working directory if necessary!
setwd('~/git/IBDpred')

require(phyloseq)
source('R/utils_to_plot.r')

# Arguments
# ===============
cohort = 'gevers' #gevers #morgan
ibd = 'cd' #uc #cd
fs = 'k40' #fcbf #k5 #k10 #k20 #k40 #deg #ldm

# Load train data
# ===============
train = readRDS('01_Preprocessing/data/train__IBD.rds')

# Load validation data
# ===============
if (cohort == 'morgan') {
  otu = read.delim2('extdata/morgan/refseq/otutable_refseq.txt', header = T, sep = '\t')
  if (ibd == 'cd') {
    clinPath = 'extdata/morgan/refseq/healthy vs cd - stool/task-healthy-cd.txt'
  } else if (ibd == 'uc'){
    clinPath = 'extdata/morgan/refseq/helthy vs uc - stool/task-healthy-uc.txt'
  }
  clin = read.delim2(clinPath, header = T, sep = '\t')
  rownames(clin) = paste0('X', clin$X.SampleID)
  clin$IBD = ifelse(clin$Var == 'Healthy', 'NO', 'YES')
  tax = readRDS('extdata/morgan/refseq/TaxonomyTable.rds')
} else if (cohort == 'gevers'){
  otu = read.delim2('extdata/gevers/refseq/otutable.txt', header = T, sep = '\t')
  clin = readRDS('extdata/gevers/task-ileum-rectum.rds')
  rownames(clin) = clin$X.SampleID
  clin$IBD = ifelse(clin$Var == 'no', 'NO', 'YES')
  tax = readRDS('extdata/gevers/refseq/TaxonomyTable.rds')
}

TAXA = tax_table(as.matrix(tax))
otu = subset(otu, select = -c(X.OTU.ID))
rownames(otu) = rownames(tax)
otu = as.matrix(otu)
OTU = otu_table(otu, taxa_are_rows = T)

CLIN = sample_data(clin)

validation = phyloseq(OTU, CLIN, TAXA)
validation = tax_glom(validation, 'Rank6')

# Intersect train and validation features
# ===============
# Get taxonomy from train and test
taxTrain = train@tax_table@.Data
taxVal = validation@tax_table@.Data

trainVars = readRDS('03_Validation/data/trainVariables.rds')
if (fs == 'fcbf') {
  signature = convertIDs(gsub('X', '', trainVars$fcbf), train)
} else if (fs == 'k5'){
  signature = convertIDs(gsub('X', '', trainVars$k5), train)
} else if (fs == 'k10'){
  signature = convertIDs(gsub('X', '', trainVars$k10), train)
} else if (fs == 'k20'){
  signature = convertIDs(gsub('X', '', trainVars$k20), train)
} else if (fs == 'k40'){
  signature = convertIDs(gsub('X', '', trainVars$k40), train)
} else if (fs == 'deg'){
  signature = convertIDs(gsub('X', '', trainVars$deg), train)
} else if(fs == 'ldm'){
  signature = convertIDs(gsub('X', '', trainVars$ldm), train)
}
otus.int = intersect(signature, taxVal[,'Rank6'])

# Extract IDs in train and validation data
# ===============
taxTrain.int = taxTrain[match(otus.int, taxTrain[,'Rank6']), ]
taxVal.int = taxVal[match(otus.int, taxVal[,'Rank6']), ]

ids.train = rownames(taxTrain.int)
ids.val = rownames(taxVal.int)

# Prune taxa in train and validation
# ===============
retrain = prune_taxa(ids.train, train)
validation = prune_taxa(ids.val, validation)

print(retrain); print(validation)

# Create data.frame to retrain the model
# ================
# In train
tt = cbind.data.frame(t(otu_table(retrain)), target = get_variable(retrain, 'IBD'))
names(tt) = make.names(names(tt))

# In validation
vv = cbind.data.frame(t(otu_table(validation)), target = get_variable(validation, 'IBD'))
names(vv) = make.names(names(vv))

# Normalization (log2)
# ===================
tt = cbind.data.frame(apply(subset(tt, select = -c(target)), 2, function(x) log2(x + 1)), target = tt$target)
vv = cbind.data.frame(apply(subset(vv, select = -c(target)), 2, function(x) log2(x + 1)), target = vv$target)

# Change names for genera names
# ===================
names(tt) = c(taxTrain.int[, 'Rank6'], 'target')
names(vv) = c(taxVal.int[, 'Rank6'], 'target')

# Checking!!
# ==================
dim(tt) ; dim(vv)
stopifnot(names(tt) == names(vv))

# Save data processed
# ==================
outDir = '03_Validation/data/'
saveRDS(tt, file = paste0(outDir, 'retrain_IBD_', cohort, '_', ibd, '_', fs))
saveRDS(vv, file = paste0(outDir, 'validation_IBD_', cohort, '_', ibd, '_', fs))

