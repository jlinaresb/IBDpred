
# Create phyloseq to Morgan and Gevers data
# ==========
require(phyloseq)
setwd('~/git/IBDpred/')

# Morgan
# ===
otu = read.delim2('extdata/morgan/refseq/otutable_refseq.txt', header = T, sep = '\t')
clinPath = 'extdata/morgan/refseq/helthy vs uc - stool/task-healthy-uc.txt'
clin = read.delim2(clinPath, header = T, sep = '\t')
rownames(clin) = paste0('X', clin$X.SampleID)
clin$IBD = ifelse(clin$Var == 'Healthy', 'NO', 'YES')
tax = readRDS('extdata/morgan/refseq/TaxonomyTable.rds')

TAXA = tax_table(as.matrix(tax))
otu = subset(otu, select = -c(X.OTU.ID))
rownames(otu) = rownames(tax)
otu = as.matrix(otu)
OTU = otu_table(otu, taxa_are_rows = T)

CLIN = sample_data(clin)

validation = phyloseq(OTU, CLIN, TAXA)
validation = tax_glom(validation, 'Rank6')

morgan = validation
saveRDS(morgan, file = 'Review/data/morgan.rds')


# Gevers
# ====
otu = read.delim2('extdata/gevers/refseq/otutable.txt', header = T, sep = '\t')
clin = readRDS('extdata/gevers/task-ileum-rectum.rds')
rownames(clin) = clin$X.SampleID
clin$IBD = ifelse(clin$Var == 'no', 'NO', 'YES')
tax = readRDS('extdata/gevers/refseq/TaxonomyTable.rds')

TAXA = tax_table(as.matrix(tax))
otu = subset(otu, select = -c(X.OTU.ID))
rownames(otu) = rownames(tax)
otu = as.matrix(otu)
OTU = otu_table(otu, taxa_are_rows = T)

CLIN = sample_data(clin)

validation = phyloseq(OTU, CLIN, TAXA)
validation = tax_glom(validation, 'Rank6')

gevers = validation
saveRDS(gevers, file = 'Review/data/gevers.rds')

