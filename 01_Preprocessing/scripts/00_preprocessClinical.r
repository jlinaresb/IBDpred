# American Gut Project Data Preprocessed
## NOTE!! Change working directory if necessary!
setwd('~/git/IBDpred')

## Required packages
require(dplyr)
require(fastDummies)
source('R/utils.r')

## Arguments
inputDir = 'extdata/'
clinFile = 'ag-cleaned.txt'
outDir = '01_Preprocessing/data/'

## Clinical Data
clinical = utils::read.delim(paste0(inputDir, clinFile), stringsAsFactors = F, sep = "\t", header = T, row.names = 1)
clinical$AGE_CORRECTED = as.vector(as.numeric(clinical$AGE_CORRECTED))
Clin = clinical[which(clinical$HMP_SITE == 'FECAL'), ]
Clin = clinical[which(!is.na(clinical$AGE_CORRECTED)),]
Clin$ages = cut(Clin$AGE_CORRECTED, breaks=c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 110),
                labels = c('<10', '10s', '20s', '30s', '40s', '50s', '60s', '70s', '80s', '90+'),
                right = FALSE)

SampleData = Clin

# Save ID patients in a variable
SampleData$IDs = rownames(SampleData)

# Delete VIOSCREEN variables
VSFeats = grep('VIOSCREEN', colnames(SampleData))
clinical = SampleData[, - VSFeats]

# Select target variables
feats = c('IDs','COUNTRY_OF_BIRTH', 'LIVER_DISEASE', 'SUBSET_HEALTHY', 'SUBSET_BMI', 'AGE_YEARS', 'EXERCISE_FREQUENCY', 
          'DOG', 'FERMENTED_PLANT_FREQUENCY', 'SUGARY_SEWWTS_FREQUENCY', 'ANONYMIZED_NAME',
          'HOMECOOKED_MEALS_FREQUENCY', 'ECONOMIC_REGION', 'NAIL_BITER', 'ADD_ADHD', 'IBD', 'IBS',
          'DIET_TYPE', 'CANCER', 'DIABETES', 'CDIFF', 'MENTAL_ILLNESS_TYPE_BIPOLAR_DISORDER', 'HEIGHT_CM', 'CHICKENPOX',
          'FUNGAL_OVERGROWTH', 'GEO_LOC_NAME', 'HOST_SUBJECT_ID', 'CARDIOVASCULAR_DISEASE', 'FED_AS_INFANT', 
          'MENTAL_ILLNESS', 'COUNTRY', 'ACID_REFLUX', 'PREGNANT', 'BMI', 'TEETHBRUSHING_FREQUENCY', 'DEPRESSION_BIPOLAR_SCHIZOPHRENIA',
          'SEX', 'ASD', 'LUNG_DISEASE', 'ALCOHOL_FREQUENCY', 'SUBSET_ANTIBIOTIC_HISTORY', 'SMOKING_FREQUENCY', 'LEVEL_OF_EDUCATION',
          'MENTAL_ILLNESS_TYPE_DEPRESSION', 'KIDNEY_DISEASE', 'COLLECTION_MONTH', 'THYROID', 'GLUTEN', 'WEIGHT_KG', 'ALLERGIC_TO_PEANUTS',
          'PROBIOTIC_FREQUENCY', 'CANCER_TREATMENT', 'MENTAL_ILLNESS_TYPE_SCHIZOPHRENIA', 'IBD_DIAGNOSIS', 'SLEEP_DURATION', 'LACTOSE',
          'BOWEL_MOVEMENT_QUALITY', 'RACE', 'ALZHEIMERS', 'ANTIBIOTIC_HISTORY', 'HIGH_FAT_RED_MEAT_FREQUENCY', 'APPENDIX_REMOVED', 
          'DIABETES_TYPE', 'BMI_CAT', 'DRINKING_WATER_SOURCE', 'RED_MEAT_FREQUENCY', 'AUTOIMMUNE', 'BIRTH_YEAR', 'ALCOHOL_CONSUMPTION',
          'VEGETABLE_FREQUENCY', 'HMP_SITE', 'ages', 'AGE_CORRECTED', 'COLLECTION_TIMESTAMP')

clinical = clinical[, intersect(feats, names(clinical))]

numF = c('AGE_YEARS', 'HEIGHT_CM', 'BMI', 'WEIGHT_KG', 'BIRTH_YEAR', 'AGE_CORRECTED')
logicalF = c('SUBSET_HEALTHY', 'SUBSET_BMI', 'DOG', 'NAIL_BITER', 'SUBSET_ANTIBIOTIC_HISTORY', 'ALLERGIC_TO_PEANUTS', 'LACTOSE', 'ALCOHOL_CONSUMPTION', 'APPENDIX_REMOVED')
charF = c('IDs', 'ANONYMIZED_NAME', 'HOST_SUBJECT_ID', 'COLLECTION_TIMESTAMP')
dateF = c('COLLECTION_TIMESTAMP')
factorF = setdiff(names(clinical), c(numF, logicalF, charF, dateF))

num = as.data.frame(apply(clinical[, intersect(numF, names(clinical))], 2, function(x) as.numeric(as.character(x))))
logical = as.data.frame(apply(clinical[, intersect(logicalF, names(clinical))], 2, function(x) ifelse(x == 'true' | x == 'True' | x == 'Yes', TRUE, 
                                                                                                      ifelse(x == 'false' | x == 'False' | x == 'No', FALSE, NA))))
char = clinical[,charF]
factor = as.data.frame(lapply(clinical[, intersect(factorF, names(clinical))], as.factor))
date = strptime(clinical[, dateF], '%m/%d/%Y %H:%M')

cc = cbind.data.frame(num, logical, char, factor)
cc$COLLECTION_TIMESTAMP = strptime(cc$COLLECTION_TIMESTAMP, '%m/%d/%Y %H:%M')

fecal = cc[which(cc$HMP_SITE == 'FECAL'), ]
res = selectDupPatients(fecal)


if (dir.exists(outDir) == FALSE) {
  dir.create(outDir)
  message('Creating outDir directory ...')
}

saveRDS(res, file = '01_Preprocessing/data/clinical_SelectFeatures_unique_patients_fecal.rds')