# Merge data Morgan and Gevers
# ===
require(phyloseq)
setwd('~/git/IBDpred/')

morgan = readRDS('Review/data/morgan.rds')
gevers = readRDS('Review/data/gevers.rds')

common = intersect(rownames(otu_table(morgan)), rownames(otu_table(gevers)))

morgan = prune_taxa(common, morgan)
gevers = prune_taxa(common, gevers)

merge = merge_phyloseq(morgan, gevers)
merge.tax = tax_table(merge)@.Data



# Train data from AGP
# ==
agp = readRDS('01_Preprocessing/data/train__IBD.rds')
agp.tax = tax_table(agp)@.Data

common.with.train = intersect(merge.tax[, 'Rank6'], agp.tax[, 'Rank6'])

merge.tax = merge.tax[match(common.with.train, merge.tax[,'Rank6']),]
agp.tax = agp.tax[match(common.with.train, agp.tax[,'Rank6']),]

merge.id = rownames(merge.tax)
agp.id = rownames(agp.tax)

merge = prune_taxa(merge.id, merge)
agp = prune_taxa(agp.id, agp)

saveRDS(merge, file = 'Review/data/merge.rds')
saveRDS(agp, file = 'Review/data/agp.rds')


