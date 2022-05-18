# Machine Learning based microbiome signature to predict Inflammatory Bowel Disease subtypes

This repository includes code based on the following publication:

*Liñares Blanco J., Fernandez-Lozano C., Seoane, JA., López-Campos, G. (2022) Machine Learning based microbiome signature to predict Inflammatory Bowel Disease subtypes. Front. Microbiology 

DOI: [10.3389/fmicb.2022.872671](https://doi.org/10.3389/fmicb.2022.872671)


## Abstract

Inflammatory bowel disease (IBD) is a chronic disease with unknown pathophysiological mechanisms. There is evidence of role of microorganims in this disease development. Because of the open access to multiple omics data, it is possible to develop predictive models that are able to prognosticate the course and development of the disease. The interpretability of these models, and the study of the variables used, allows the identification of biological aspects of great importance in the development of the disease. 
In this work we obtained a metagenomic signature with predictive capacity according to the presence of IBD in faecal samples. Different Machine Learning models were trained, obtaining high performance measures. The predictive capacity of the signature was validated in two external cohorts. Specifically, in one cohort containing samples from Ulcerative Colitis and one from Crohn's Disease, the two major subtypes of IBD. The results obtained in this validation (AUC 0.74 and AUC=0.76, respectively) show that our signature presents a generalisation capacity in both subtypes. The study of the variables within the model, and a correlation study based on text mining, identified different genera that play an important and common role in the development of these two subtypes.



## Prerequisites:

Before run the code, make sure you have the following packages we've used:

```{r}
# CRAN packages
install.packages(c('dplyr', 'tidyverse', 'fastDummies', 'h2o', 'FCBF', 'mlr', 'ggplot2', 'ggpubr', 'grid', 'ggupset', 'ggalt', 'reshape2', 'forestplot', 'viridis', 'ROCR'))

# Bioconductor packages
BiocManager::install(c('phyloseq'))

```


## Project workflow

First, you must download AGP from this FTP server: ftp://ftp.microbio.me/. Navigate to ```/latest``` directory and download ```03_otus.zip``` file. In addition, extract ```ag-cleanead.txt``` file from ```04-meta.zip``` folder. To run this script, storage these data in ```/extdata``` folder.

### 01. Data preprocessing

Run by order the following scripts:

1. [00_preprocessClinical.r](https://github.com/jlinaresb/IBDpred/blob/master/01_Preprocessing/scripts/00_preprocessClinical.r)
This script have as an input the clinical data downloaded from AGP and make several processing steps in order to select and format target variables. 

2. [01_create_phyloseq.r](https://github.com/jlinaresb/IBDpred/blob/master/01_Preprocessing/scripts/01_create_phyloseq.r)
With ```.biom```, ```.tree``` and ```sampleData``` as a input, this script generates ```phyloseq``` object to manage the microbiome data. For more info to ```phyloseq``` objects, consult this [link](https://joey711.github.io/phyloseq/).

3. [02_preprocess_phyloseq.r](https://github.com/jlinaresb/IBDpred/blob/master/01_Preprocessing/scripts/02_preprocess_phyloseq.r)
After ```phyloseq``` object generated, we carried out several preprocessing steps before models trainning. First, OTUs were glomerated by genus, we performed outliers analysis through Isolation Forest algorithm using [h2o](https://cran.r-project.org/web/packages/h2o/index.html) r package before to label and balance input data. The second step was split data into train and test subset (85% train and 15% test). Only with train data, we carried out feature selection techniques. Finally, this script will generate ```.rds``` that will be the input of both Machine Learning algorithms.


### 02. Models trainnig

Run the [sh-generator.r](https://github.com/jlinaresb/IBDpred/blob/master/02_Training/scripts/sh-generator.r) script. This script generates one ```.sh``` file for each data:algorithm combination. In this case, we generated seven input datasets corresponding to each Feature Selection used (K5, K10, K20, K40, FCBF, DEG and LDM) and ran two different ML algorithms (glmnet and RF). In total, fourteen '''.sh''' files were generated. Therefore, fourteen ```.rds``` files are stored in ```02_Training/results/```.

Note that these scripts are adapted to run in HPC environment, specifically in [CESGA](https://www.cesga.es/en/home-2/). You must adapt the [config-file.r](https://github.com/jlinaresb/IBDpred/blob/master/02_Training/scripts/config-file.r) and, if necessary, the [sh-generation.r](https://github.com/jlinaresb/IBDpred/blob/master/02_Training/scripts/sh-generator.r) to your HPC environment.


### 03. Models validation
Firts, see [datasets.md](https://github.com/jlinaresb/IBDpred/blob/master/datasets.md) file for instructions to download external validation datasets. 

As we commented in the manuscript, we were not able to obtain the same features in validation data. As a result we needed to retrain models only with available features. 

First, we need to extract features of each input dataset and save them in a ```list``` object. [get_train_variables.r](https://github.com/jlinaresb/IBDpred/blob/master/03_Validation/scripts/get_train_variables.r) script make this step.

After that, you can run the [external_validation.r](https://github.com/jlinaresb/IBDpred/blob/master/03_Validation/scripts/external_validation.r) script. Change the variables in argument sections to obtain new train and validation cohorts.

Now, [sh-generator.r](https://github.com/jlinaresb/IBDpred/blob/master/02_Training/scripts/sh-generator.r) can be run again after change arguments in [config-file.r](https://github.com/jlinaresb/IBDpred/blob/master/02_Training/scripts/config-file.r). To follow the analysis, we recomend save ```bmr``` objects in ```03_Validation/results/retrain/```


## Figures
[Figures](https://github.com/jlinaresb/IBDpred/tree/master/Figures) folder contains scripts to reproduce each figures of the paper. Note that for aesthetics reasons, the figures were edited using illustrator.


## Review
All modifications after peer-review process were stored in [Review](https://github.com/jlinaresb/IBDpred/tree/master/Review folder. 


## Citation:
```{tex}
@ARTICLE{10.3389/fmicb.2022.872671,
  
	AUTHOR={Liñares-Blanco, Jose and Fernandez-Lozano, Carlos and Seoane, Jose A. and López-Campos, Guillermo},   
		 
	TITLE={Machine Learning Based Microbiome Signature to Predict Inflammatory Bowel Disease Subtypes},      
		
	JOURNAL={Frontiers in Microbiology},      
		
	VOLUME={13},      
		
	YEAR={2022},      
		  
	URL={https://www.frontiersin.org/article/10.3389/fmicb.2022.872671},       
		
	DOI={10.3389/fmicb.2022.872671},      
		
	ISSN={1664-302X},   
	   
	ABSTRACT={Inflammatory bowel disease (IBD) is a chronic disease with unknown pathophysiological mechanisms. There is evidence of the role of microorganims in this disease development. Thanks to the open access to multiple omics data, it is possible to develop predictive models that are able to prognosticate the course and development of the disease. The interpretability of these models, and the study of the variables used, allows the identification of biological aspects of great importance in the development of the disease. In this work we generated a metagenomic signature with predictive capacity to identify IBD from fecal samples. Different Machine Learning models were trained, obtaining high performance measures. The predictive capacity of the identified signature was validated in two external cohorts. More precisely a cohort containing samples from patients suffering Ulcerative Colitis and another from patients suffering Crohn's Disease, the two major subtypes of IBD. The results obtained in this validation (AUC 0.74 and AUC = 0.76, respectively) show that our signature presents a generalization capacity in both subtypes. The study of the variables within the model, and a correlation study based on text mining, identified different genera that play an important and common role in the development of these two subtypes.}
}
```

## Questions?
If you have any questions, please feel free to contact (jose.linares@genyo.es).