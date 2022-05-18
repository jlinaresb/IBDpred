glom_outliers = function(filePath, taxLevel){
  
  require(phyloseq)
  source('~/git/Metagenomics/isolationForest.r')
  
  # Load phyloseq object
  print('Reading phyloseq object...')
  pseq = readRDS(filePath)
  
  # Glomerate by taxonomy
  print(paste0('Glomerating OTUs by ', taxLevel, ' ...'))
  # pseqGlom = tax_glom(pseq, taxLevel)
  
  # Outlier analysis (Isolation Forest)
  res = isoForest(pseq)
  print(res)
  
  return(res)
}



labeling = function(physeq, target){
  
  require(phyloseq)
  
  print(paste0('Labeling data by ', target))
  
  if (target == 'LUNG_DISEASE') {
    physeq = subset_samples(physeq, LUNG_DISEASE == 'Diagnosed by a medical professional (doctor, physician assistant)' |
                              LUNG_DISEASE == 'I do not have this condition')
    Y = sample_data(physeq)$LUNG_DISEASE
    Y = ifelse(Y == 'I do not have this condition', 'NO', 'YES')
    sample_data(physeq)$LUNG_DISEASE = Y
    # Balance Data
    major = rownames(sample_data(subset_samples(physeq, LUNG_DISEASE == 'NO')))
    minor = rownames(sample_data(subset_samples(physeq, LUNG_DISEASE == 'YES')))
    set.seed(5555)
    subSampling = sample(x = major, size = length(minor))
    physeq = prune_samples(c(subSampling, minor), physeq)
    
  }else if (target == 'IBD'){
    physeq = subset_samples(physeq, IBD == 'Diagnosed by a medical professional (doctor, physician assistant)' |
                              IBD == 'I do not have this condition')
    Y = sample_data(physeq)$IBD
    Y = ifelse(Y == 'I do not have this condition', 'NO', 'YES')
    sample_data(physeq)$IBD = Y
    # Balance Data
    major = rownames(sample_data(subset_samples(physeq, IBD == 'NO')))
    minor = rownames(sample_data(subset_samples(physeq, IBD == 'YES')))
    set.seed(5555)
    subSampling = sample(x = major, size = length(minor))
    physeq = prune_samples(c(subSampling, minor), physeq)
    
  }else if (target == 'IBS'){
    physeq = subset_samples(physeq, IBS == 'Diagnosed by a medical professional (doctor, physician assistant)' |
                                 IBS == 'I do not have this condition')
    Y = sample_data(physeq)$IBS
    Y = ifelse(Y == 'I do not have this condition', 'NO', 'YES')
    sample_data(physeq)$IBS = Y
    # Balance Data
    major = rownames(sample_data(subset_samples(physeq, IBS == 'NO')))
    minor = rownames(sample_data(subset_samples(physeq, IBS == 'YES')))
    set.seed(5555)
    subSampling = sample(x = major, size = length(minor))
    physeq = prune_samples(c(subSampling, minor), physeq)
    
  }else if (target == 'ASD'){
    physeq = subset_samples(physeq, ASD == 'Diagnosed by a medical professional (doctor, physician assistant)' |
                              ASD == 'I do not have this condition')
    Y = sample_data(physeq)$ASD
    Y = ifelse(Y == 'I do not have this condition', 'NO', 'YES')
    sample_data(physeq)$ASD = Y
    # Balance Data
    major = rownames(sample_data(subset_samples(physeq, ASD == 'NO')))
    minor = rownames(sample_data(subset_samples(physeq, ASD == 'YES')))
    set.seed(5555)
    subSampling = sample(x = major, size = length(minor))
    physeq = prune_samples(c(subSampling, minor), physeq)
    
  }else if(target == 'DIET_TYPE'){
    physeq = subset_samples(physeq, DIET_TYPE == 'Omnivore' | DIET_TYPE == 'Vegetarian' | DIET_TYPE == 'Vegan')
    Y = sample_data(physeq)$DIET_TYPE 
    Y = ifelse(Y == 'Omnivore', 'Omnivore', 'Vegetarian/Vegan')
    sample_data(physeq)$DIET_TYPE = Y
    # Balance Data
    major = rownames(sample_data(subset_samples(physeq, DIET_TYPE == 'Omnivore')))
    minor = rownames(sample_data(subset_samples(physeq, DIET_TYPE == 'Vegetarian/Vegan')))
    set.seed(5555)
    subSampling = sample(x = major, size = length(minor))
    physeq = prune_samples(c(subSampling, minor), physeq)
    
  }else if (target == 'COUNTRY'){
    physeq = subset_samples(physeq, COUNTRY == 'USA' | COUNTRY == 'United Kingdom')
    # Balance Data
    major = rownames(sample_data(subset_samples(physeq, COUNTRY == 'USA')))
    minor = rownames(sample_data(subset_samples(physeq, COUNTRY == 'United Kingdom')))
    set.seed(5555)
    subSampling = sample(x = major, size = length(minor))
    physeq = prune_samples(c(subSampling, minor), physeq)
    
  }else if (target == 'ANTIBIOTIC_HISTORY'){
    physeq = subset_samples(physeq, ANTIBIOTIC_HISTORY == '6 months' |
                              ANTIBIOTIC_HISTORY == 'Week' |
                              ANTIBIOTIC_HISTORY == 'Month' |
                              ANTIBIOTIC_HISTORY == 'I have not taken antibiotics in the past year.')
    Y = sample_data(physeq)$ANTIBIOTIC_HISTORY
    Y = ifelse(Y == 'I have not taken antibiotics in the past year.', 'noAB-1y', 'AB-6m')
    sample_data(physeq)$ANTIBIOTIC_HISTORY = Y
    # Balance Data
    major = rownames(sample_data(subset_samples(physeq, ANTIBIOTIC_HISTORY == 'noAB-1y')))
    minor = rownames(sample_data(subset_samples(physeq, ANTIBIOTIC_HISTORY == 'AB-6m')))
    set.seed(5555)
    subSampling = sample(x = major, size = length(minor))
    physeq = prune_samples(c(subSampling, minor), physeq)
    
  }
  
  print('Balancing data ...')
  print(table(get_variable(physeq, target)))
  
  return(physeq)
}



split.physeq = function(physeq, ratio = 0.85){
  
  require(phyloseq)
  
  smp_size = floor(ratio * nsamples(physeq))
  
  set.seed(123)
  train_ind = sample(seq_len(nsamples(physeq)), size = smp_size)
  trainPats = sample_names(physeq)[train_ind]
  testPats = sample_names(physeq)[-train_ind]
  
  train = subset_samples(physeq, sample_names(physeq) %in% trainPats)
  test = subset_samples(physeq, sample_names(physeq) %in% testPats)
  
  res = list(train = train, test = test)
  
  print(paste0('Data splitted into train and test with ', ratio, ' samples in train'))
  
  return(res)
  
}



deg = function(physeq, target){
  
  require(phyloseq)
  source('~/git/Metagenomics/utils.r')
  
  #ps_relAb = transform_sample_counts(physeq, function(x) x / sum(x))
  
  ps.nzv = nzv(physeq, threshold = 50)
  dge = phyloseq_to_edgeR(ps.nzv, group = target)
  
  et = exactTest(dge)
  tt = topTags(et, n=nrow(dge$table), adjust.method="BH", sort.by="PValue")
  res = tt@.Data[[1]]
  alpha = 0.001
  sigtab = res[(res$FDR < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq)[rownames(sigtab), ], "matrix"))
  
  res = prune_taxa(rownames(sigtab), ps.nzv)
  
  print('Differential OTUs expression done!')
  
  return(res)
  
}


abRel = function(pseq){
  
  require(phyloseq)
  
  print('Transforming counts in relative abundance...')
  pseq2 = transform_sample_counts(pseq, function(x) x / sum(x))
  
  return(pseq2)
}



nzv = function(physeq, threshold){
  
  require(phyloseq)
  
  keepOTUs = names(which(rowSums(otu_table(physeq) == 0) / ncol(otu_table(physeq)) * 100 < threshold))
  physeq = subset_taxa(physeq, (Rank6 != 'g__'))
  pseq = prune_taxa(keepOTUs, physeq)
  
  print(paste0(ntaxa(physeq) - ntaxa(pseq), ' taxa were removed!'))
  
  return(pseq)
}



featureSelection = function(dataPath, outDir){
  
  require(FCBF)
  source('~/git/Metagenomics/fs_ml_functions.r')
  require(phyloseq)
  
  data = readRDS(dataPath)
  
  
  
  fcbf = fast.corr.FS(data, thresh = 0.0025)
  saveRDS(fcbf, paste0(outDir, 'fcbf_', ncol(fcbf) - 1, gsub('.rds', '', basename(dataPath))))
  fs.abs(data, 'kruskal.test', c(5, 15, 45, 135), outDir, gsub('.rds', '', basename(dataPath)))
  
  print(paste0('Feature Selection datasets saved them in ', outDir))
  
}



# Final Function
# ===============================
processing = function(filePath, outDir, taxLevel = 'Rank6', target = 'DIET_TYPE', threshold = 50, ratio = 0.85){
  
  pseq = glom_outliers(filePath, taxLevel)
  pseq = labeling(pseq, target)
  return(pseq)
  
  # split = split.physeq(pseq, ratio)
  # 
  # train = split$train
  # test = split$test
  # 
  # diffExp = deg(train, target)
  # 
  # physeq_nzv = nzv(diffExp, threshold)
  # 
  # featureSelection(physeq_nzv, outDir)
  
  print('Analysis done!')
}










