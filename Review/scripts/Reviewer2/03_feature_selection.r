require(phyloseq)
setwd('~/git/IBDpred/')
source('R/fs_ml_functions.r')

# Load data
# ===
merge = readRDS('Review/data/merge.rds')
agp = readRDS('Review/data/agp.rds')

merge.genus = unname(tax_table(merge)@.Data[,'Rank6'])
agp.genus = unname(tax_table(agp)@.Data[,'Rank6'])

counts.merge = t(otu_table(merge)@.Data)
counts.merge = apply(counts.merge, 2, function(x) log2(x + 1))
colnames(counts.merge) = merge.genus
merge = cbind.data.frame(
  counts.merge, 
  target = get_variable(merge, 'IBD')
)
saveRDS(merge, file = '~/git/IBDpred/Review/data/merge_df.rds')

counts.agp = t(otu_table(agp)@.Data)
counts.agp = apply(counts.agp, 2, function(x) log2(x + 1))
colnames(counts.agp) = agp.genus
agp = cbind.data.frame(
  counts.agp, 
  target = get_variable(agp, 'IBD')
)
saveRDS(agp, file = '~/git/IBDpred/Review/data/agp_df.rds')

# Feature Selection
# ====
outDir = '~/git/IBDpred/Review/data/fs/'
fcbf = fast.cor.FS(merge, thresh = 0.00025)
saveRDS(fcbf, file = paste0(outDir, 'fcbf_merge'))
fs.abs(merge, 'kruskal.test', c(5, 10, 20, 40), outDir, 'kruskal_merge')
