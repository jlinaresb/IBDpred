## --- Create Phyloseq object --- ##

## NOTE!! Change working directory if necessary!
setwd('~/git/IBDpred')

# Required packages
require(phyloseq)

# Defining biom, sample and tree files (100nt) --> trimmed sequencing
# NOTE!!! Sample data has been preprocessed to add 'ages' variable and subset clinical features
biom_file = 'extdata/03-otus/100nt/gg-13_8-97-percent/otu_table.biom'
tree_file = 'extdata/03-otus/100nt/gg-13_8-97-percent/97_otus.tree')
sampleDir = '01_Preprocessing/data/clinical_SelectFeatures_unique_patients_fecal.rds'

# Reading the three files into phyloseq compatible objects
print('Importing OTU table ... ')
otu = import_biom(biom_file)
print('Reading tree file ...')
tree = read_tree(tree_file)
print('Reading sample data ...')

sample = readRDS(sampleDir)
sample_cohort = sample_data(sample)

phyloseq = merge_phyloseq(otu, tree, sample_cohort)
saveRDS(phyloseq, file = '01_Preprocessing/data/phyloseq.rds')
