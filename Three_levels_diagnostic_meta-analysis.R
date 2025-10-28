# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - Three-Level Model for Diagnostic Meta-Analysis
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load required libraries
library(readxl)
library(metafor)
library(ggplot2)

# Set working directory
setwd("data_path")
# Load data
data <- read_excel("data.xlsx","clean_data")

# Effect size calculation
log_dor <- escalc(measure = "OR", ai = tp, bi = fn, ci = fp, di = tn, n1i = PD, n2i = HC, data = data, add = 1/2, to = "all")

#Three level meta-analysis for DOR
three_level_dor <- rma.mv(yi, vi, random = ~ 1 | Study_num/Table_num, data = log_dor, slab = `Study`)
log_DOR_predicted <- predict(three_level_dor, digits=2)
DOR_predicted <- predict(three_level_dor, transf=exp, digits=2)
print(DOR_predicted, digits=2)
I_sqr <- var.comp(three_level_dor) # To run var.comp function, one need to run the below code in Console [https://raw.githubusercontent.com/MathiasHarrer/dmetar/master/R/mlm.variance.distribution.R]
print(I_sqr)
summary(three_level_dor)

#Moderator analysis
log_dor$Feature_Selection <- relevel(factor(log_dor$Feature_Selection), ref = "Only GTA") 
log_dor$Algorithm <- relevel(factor(log_dor$Algorithm), ref = "ROC Analysis")
log_dor$Imaging_Technology <- relevel(factor(log_dor$Imaging_Technology), ref = "Structural Imaging")
log_dor$Threshold <- relevel(factor(log_dor$Threshold), ref = "proportion")
control_feature_log_dor <- log_dor[log_dor$Feature_Selection == "Only GTA", ]

mods_analysis_algorithm <- rma.mv(yi, vi, mods = ~factor(`Algorithm`), random = ~ 1| `Study_num`/`Table_num`,data=log_dor,slab=`Study`)
summary(mods_analysis_algorithm)

mods_analysis_imaging_technology <- rma.mv(yi, vi, mods = ~factor(`Imaging_Technology`), random = ~ 1| `Study_num`/`Table_num`,data=control_feature_log_dor,slab=`Study`)
summary(mods_analysis_imaging_technology)

mods_analysis_feature_selection <- rma.mv(yi, vi, mods = ~factor(`Feature_Selection`), random = ~ 1| `Study_num`/`Table_num`,data=log_dor,slab=`Study`)
summary(mods_analysis_feature_selection)

mods_analysis_threshold <- rma.mv(yi, vi, mods = ~factor(`Threshold`), random = ~ 1| `Study_num`/`Table_num`,data=log_dor,slab=`Study`)
summary(mods_analysis_threshold)
