
# set working directory 
setwd(rstudioapi::selectDirectory(caption = "Please choose the folder, where the R script files were stored"))


library(pacman)
pacman::p_load(readxl,dplyr,ggplot2,tibble,metafor,googledrive,tidyr, data.table,cowplot,
               extrafont,scales,showtext,colorspace,purrr,ggbeeswarm,stringr)


source("function.R")


###############################################dMRI#####################
modality="dMRI"
source("data_processing.R")
meta.d=meta_d(modality)
############## Clustering coefficient################
metric="Cp"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm


#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger


#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim

inf=influence(unim)
plot(inf)

# Remove influential effect size (Nigro2016)
data=data[-20,]

#re-run multiple meta after remove influence study
mlm.inf <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm.inf


###meta-regression
data=meta.d[[metric]]

##Continuous variable
continuous_vars = c("PD.Age.years", "PD.Male.prop", "PD.Education.years", "OnsetAge", 
                    "Duration.years", "UPDRS_III", "H_Y", "MMSE","LEDD","Nodes")

continuous_results = cont_res(continuous_vars)

##Categorical variable
categorical_vars = c("Drug_code", "Template_AAL","Template_cerebellum","Tractography_code","Edges.weighting_code", "Threshold_code")
data=process_categorical_vars(categorical_vars, data)

categorical_results = cate_res(categorical_vars,data)

dmri=list()
dmri.d=c(list(
  mlm=mlm,
  egger=egger,
  mlm.inf=mlm.inf
),continuous_results,categorical_results)

dmri[[metric]]=dmri.d
rm(list = setdiff(ls(), c("meta.d", "dmri", lsf.str())))


############## Local efficiency###############
metric="Eloc"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm

#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim

inf=influence(unim)
plot(inf)

##Continuous variable
continuous_vars = c("PD.Age.years", "PD.Male.prop", "PD.Education.years", "OnsetAge", 
                    "Duration.years", "UPDRS_III", "H_Y", "MMSE","LEDD","Nodes")

continuous_results = cont_res(continuous_vars)

##Categorical variable
categorical_vars = c("Drug_code", "Template_AAL","Template_cerebellum","Tractography_code","Edges.weighting_code", "Threshold_code")
data=process_categorical_vars(categorical_vars, data)

categorical_results = cate_res(categorical_vars,data)

dmri.d=c(list(
  mlm=mlm,
  egger=egger
),continuous_results,categorical_results)

dmri[[metric]]=dmri.d
rm(list = setdiff(ls(), c("meta.d", "dmri", lsf.str())))



############## Characteristic path length################
metric="Lp"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm

#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim

inf=influence(unim)
plot(inf)

#remove 7 effect size(Gan et al.2024 PD-NICD)
data=data[-7,]

#re-run multiple meta after remove influence study
mlm.inf <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm.inf

###meta-regression
data=meta.d[[metric]]

##Continuous variable
continuous_vars = c("PD.Age.years", "PD.Male.prop", "PD.Education.years", "OnsetAge", 
                    "Duration.years", "UPDRS_III", "H_Y", "MMSE","LEDD","Nodes")

continuous_results = cont_res(continuous_vars)

##Categorical variable
categorical_vars = c("Drug_code", "Template_AAL","Template_cerebellum","Tractography_code","Edges.weighting_code", "Threshold_code")
data=process_categorical_vars(categorical_vars, data)

categorical_results = cate_res(categorical_vars,data)

dmri.d=c(list(
  mlm=mlm,
  egger=egger,
  mlm.inf=mlm.inf
),continuous_results,categorical_results)

dmri[[metric]]=dmri.d
rm(list = setdiff(ls(), c("meta.d", "dmri", lsf.str())))

############## Global efficiency################
metric="Eglob"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm

#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim

inf=influence(unim)
plot(inf)

##Continuous variable
continuous_vars = c("PD.Age.years", "PD.Male.prop", "PD.Education.years", "OnsetAge", 
                    "Duration.years", "UPDRS_III", "H_Y", "MMSE","LEDD","Nodes")

continuous_results = cont_res(continuous_vars)

##Categorical variable
categorical_vars = c("Drug_code","Template_AAL","Template_cerebellum","Tractography_code","Edges.weighting_code", "Threshold_code")
data=process_categorical_vars(categorical_vars, data)

categorical_results = cate_res(categorical_vars,data)

dmri.d=c(list(
  mlm=mlm,
  egger=egger
),continuous_results,categorical_results)

dmri[[metric]]=dmri.d
rm(list = setdiff(ls(), c("meta.d", "dmri", lsf.str())))

############## normalized clustering coefficient################
metric="Gamma"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm


#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim

inf=influence(unim)
plot(inf)


##Continuous variable
continuous_vars = c("PD.Age.years", "PD.Male.prop", "PD.Education.years", "OnsetAge", 
                    "Duration.years", "UPDRS_III", "H_Y", "MMSE","LEDD","Nodes")

continuous_results = cont_res(continuous_vars)

##Categorical variable
categorical_vars = c( "Drug_code", "Template_AAL","Template_cerebellum","Tractography_code","Edges.weighting_code", "Threshold_code")
data=process_categorical_vars(categorical_vars, data)

categorical_results = cate_res(categorical_vars,data)

dmri.d=c(list(
  mlm=mlm,
  egger=egger
),continuous_results,categorical_results)

dmri[[metric]]=dmri.d
rm(list = setdiff(ls(), c("meta.d", "dmri", lsf.str())))


############## normalized Characteristic path length################
metric="Lambda"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm

#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim

inf=influence(unim)
plot(inf)


##Continuous variable
continuous_vars = c("PD.Age.years", "PD.Male.prop", "PD.Education.years", "OnsetAge", 
                    "Duration.years", "UPDRS_III", "H_Y", "MMSE","LEDD","Nodes")

continuous_results = cont_res(continuous_vars)

##Categorical variable
categorical_vars = c( "Drug_code","Template_AAL","Template_cerebellum","Tractography_code","Edges.weighting_code", "Threshold_code")
data=process_categorical_vars(categorical_vars, data)

categorical_results = cate_res(categorical_vars,data)

dmri.d=c(list(
  mlm=mlm,
  egger=egger
),continuous_results,categorical_results)

dmri[[metric]]=dmri.d
rm(list = setdiff(ls(), c("meta.d", "dmri", lsf.str())))

############## small-worldness################
metric="Sigma"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm


#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim

inf=influence(unim)
plot(inf)

##Continuous variable
continuous_vars = c("PD.Age.years", "PD.Male.prop", "PD.Education.years", "OnsetAge", 
                    "Duration.years", "UPDRS_III", "H_Y", "MMSE","LEDD","Nodes")

continuous_results = cont_res(continuous_vars)

##Categorical variable
categorical_vars = c("Drug_code", "Template_AAL","Template_cerebellum","Tractography_code","Edges.weighting_code", "Threshold_code")
data=process_categorical_vars(categorical_vars, data)

categorical_results = cate_res(categorical_vars,data)


dmri.d=c(list(
  mlm=mlm,
  egger=egger
),continuous_results,categorical_results)

dmri[[metric]]=dmri.d
rm(list = setdiff(ls(), c("meta.d", "dmri", lsf.str())))


############## Transitivity################
metric="Transitivity"
data=meta.d[[metric]]

# uni meta
unim <- rma.uni(yi, vi, slab = Study, data = data)
unim

dmri.d=list(
  unim=unim
)

dmri[[metric]]=dmri.d
rm(list = setdiff(ls(), c("meta.d", "dmri", lsf.str())))


############## Modularity################
metric="Modularity"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm


#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim

inf=influence(unim)
plot(inf)

#remove 3 effect size (Vriend et al. 2018)
data=data[-3,]

#re-run multiple meta after remove influence study
mlm.inf <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm.inf

dmri.d=list(
  mlm=mlm,
  egger=egger,
  mlm.inf=mlm.inf
)

dmri[[metric]]=dmri.d
rm(list = setdiff(ls(), c("meta.d", "dmri", lsf.str())))



############## Assortativity################
metric="Assortativity"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm

dmri.d=list(
  mlm=mlm
)

dmri[[metric]]=dmri.d
rm(list = setdiff(ls(), c("meta.d", "dmri", lsf.str())))

############## BC################
metric="BC"
data=meta.d[[metric]]

# uni meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm

dmri.d=list(
  mlm=mlm
)

dmri[[metric]]=dmri.d
rm(list = setdiff(ls(), c("meta.d", "dmri", lsf.str())))

############## Strength################
metric="Strength"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm


#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim

inf=influence(unim)
plot(inf)

dmri.d=list(
  mlm=mlm,
  egger=egger
)

dmri[[metric]]=dmri.d
rm(list = setdiff(ls(), c("meta.d", "dmri", lsf.str())))



############## Degree################
metric="Degree"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm


#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim

inf=influence(unim)
plot(inf)

dmri.d=list(
  mlm=mlm,
  egger=egger
)

dmri[[metric]]=dmri.d
rm(list = setdiff(ls(), c("meta.d", "dmri", lsf.str())))


############## Density################
metric="Density"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm


#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim

inf=influence(unim)
plot(inf)

dmri.d=list(
  mlm=mlm,
  egger=egger
)

dmri[[metric]]=dmri.d
rm(list = setdiff(ls(), c("meta.d", "dmri", lsf.str())))



############## Eccentricity################
metric="Eccentricity"
data=meta.d[[metric]]

# uni meta
unim <- rma.uni(yi, vi, slab = Study, data = data)
unim

dmri.d=list(
  unim=unim
)

dmri[[metric]]=dmri.d
rm(list = setdiff(ls(), c("meta.d", "dmri", lsf.str())))



############## save data#########
results=list()
results[["dMRI"]]=dmri
###############################################fMRI####################
modality="fMRI"
source("data_processing.R")
meta.d=meta_d(modality)
############## Clustering coefficient######
metric="Cp"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm


#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim

inf=influence(unim)
plot(inf)

##Continuous variable
continuous_vars = c("PD.Age.years", "PD.Male.prop", "PD.Education.years", "OnsetAge", 
                    "Duration.years", "UPDRS_III", "H_Y", "MMSE","LEDD","Nodes")

continuous_results = cont_res(continuous_vars)

##Categorical variable
categorical_vars = c("Drug_code", "Scanning", "Template_code","Network_framework", "Threshold_code")
data=process_categorical_vars(categorical_vars, data)

categorical_results = cate_res(categorical_vars,data)

fmri=list()
fmri.d=c(list(
  mlm=mlm,
  egger=egger
),continuous_results,categorical_results)

fmri[[metric]]=fmri.d
rm(list = setdiff(ls(), c("meta.d", "fmri","results", lsf.str())))

############## Local efficiency######
metric="Eloc"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm


#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim

inf=influence(unim)
plot(inf)

##Continuous variable
continuous_vars = c("PD.Age.years", "PD.Male.prop", "PD.Education.years", "OnsetAge", 
                    "Duration.years", "UPDRS_III", "H_Y", "MMSE","LEDD","Nodes")

continuous_results = cont_res(continuous_vars)

##Categorical variable
categorical_vars = c("Drug_code", "Scanning", "Template_code","Network_framework", "Threshold_code")
data=process_categorical_vars(categorical_vars, data)

categorical_results = cate_res(categorical_vars,data)

fmri.d=c(list(
  mlm=mlm,
  egger=egger
),continuous_results,categorical_results)

fmri[[metric]]=fmri.d
rm(list = setdiff(ls(), c("meta.d", "fmri","results", lsf.str())))

############## Characteristic path length######
metric="Lp"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm


#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim

inf=influence(unim)
plot(inf)

##Continuous variable
continuous_vars = c("PD.Age.years", "PD.Male.prop", "PD.Education.years", "OnsetAge", 
                    "Duration.years", "UPDRS_III", "H_Y", "MMSE","LEDD","Nodes")

continuous_results = cont_res(continuous_vars)

##Categorical variable
categorical_vars = c("Drug_code", "Scanning", "Template_code","Network_framework", "Threshold_code")
data=process_categorical_vars(categorical_vars, data)

categorical_results = cate_res(categorical_vars,data)

fmri.d=c(list(
  mlm=mlm,
  egger=egger
),continuous_results,categorical_results)

fmri[[metric]]=fmri.d
rm(list = setdiff(ls(), c("meta.d", "fmri","results", lsf.str())))
############## Global efficiency######
metric="Eglob"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm


#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim

inf=influence(unim)
plot(inf)

##remove 22 effect size(Lopes et al. 2017 Group 4) 

data=data[-22,]

#re-run multiple meta after remove influence study
mlm.inf <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm.inf


###meta-regression
data=meta.d[[metric]]

##Continuous variable
continuous_vars = c("PD.Age.years", "PD.Male.prop", "PD.Education.years", "OnsetAge", 
                    "Duration.years", "UPDRS_III", "H_Y", "MMSE","LEDD","Nodes")

continuous_results = cont_res(continuous_vars)

##Categorical variable
categorical_vars = c("Drug_code", "Scanning", "Template_code","Network_framework", "Threshold_code")
data=process_categorical_vars(categorical_vars, data)

categorical_results = cate_res(categorical_vars,data)


fmri.d=c(list(
  mlm=mlm,
  egger=egger,
  mlm.inf=mlm.inf
),continuous_results,categorical_results)

fmri[[metric]]=fmri.d
rm(list = setdiff(ls(), c("meta.d", "fmri","results", lsf.str())))
############## normalized clustering coefficient######
metric="Gamma"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm


#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim

inf=influence(unim)
plot(inf)

##Continuous variable
continuous_vars = c("PD.Age.years", "PD.Male.prop", "PD.Education.years", "OnsetAge", 
                    "Duration.years", "UPDRS_III", "H_Y", "MMSE","LEDD","Nodes")

continuous_results = cont_res(continuous_vars)

##Categorical variable
categorical_vars = c("Drug_code", "Scanning", "Template_code","Network_framework", "Threshold_code")
data=process_categorical_vars(categorical_vars, data)

categorical_results = cate_res(categorical_vars,data)

fmri.d=c(list(
  mlm=mlm,
  egger=egger
),continuous_results,categorical_results)

fmri[[metric]]=fmri.d
rm(list = setdiff(ls(), c("meta.d", "fmri","results", lsf.str())))
############## normalized Characteristic path length######
metric="Lambda"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm


#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim

inf=influence(unim)
plot(inf)

##remove 8 effect size(	Shang et al. 2023b)

data=data[-8,]

#re-run multiple meta after remove influence study
mlm.inf <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm.inf


###meta-regression
data=meta.d[[metric]]

##Continuous variable
continuous_vars = c("PD.Age.years", "PD.Male.prop", "PD.Education.years", "OnsetAge", 
                    "Duration.years", "UPDRS_III", "H_Y", "MMSE","LEDD","Nodes")

continuous_results = cont_res(continuous_vars)

##Categorical variable
categorical_vars = c("Drug_code", "Scanning", "Template_code","Network_framework", "Threshold_code")
data=process_categorical_vars(categorical_vars, data)

categorical_results = cate_res(categorical_vars,data)

fmri.d=c(list(
  mlm=mlm,
  egger=egger,
  mlm.inf=mlm.inf
),continuous_results,categorical_results)

fmri[[metric]]=fmri.d
rm(list = setdiff(ls(), c("meta.d", "fmri","results", lsf.str())))
############## small-worldness######
metric="Sigma"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm


#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim

inf=influence(unim)
plot(inf)

##remove 9 effect size(Siva et al. 2024PD-ODP)
data=data[-9,]

#re-run multiple meta after remove influence study
mlm.inf <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm.inf


###meta-regression
data=meta.d[[metric]]

##Continuous variable
continuous_vars = c("PD.Age.years", "PD.Male.prop", "PD.Education.years", "OnsetAge", 
                    "Duration.years", "UPDRS_III", "H_Y", "MMSE","LEDD","Nodes")

continuous_results = cont_res(continuous_vars)

##Categorical variable
categorical_vars = c("Drug_code", "Scanning", "Template_code","Network_framework", "Threshold_code")
data=process_categorical_vars(categorical_vars, data)

categorical_results = cate_res(categorical_vars,data)


fmri.d=c(list(
  mlm=mlm,
  egger=egger,
  mlm.inf=mlm.inf
),continuous_results,categorical_results)

fmri[[metric]]=fmri.d
rm(list = setdiff(ls(), c("meta.d", "fmri","results", lsf.str())))
############## transitivity######
metric="Transitivity"
data=meta.d[[metric]]

# meta
unim <- rma.uni(yi, vi, slab = Study, data = data, method = "REML")
unim


fmri.d=list(
  unim=unim
)

fmri[[metric]]=fmri.d
rm(list = setdiff(ls(), c("meta.d", "fmri","results", lsf.str())))

############## Modularity######
metric="Modularity"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm


#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim

inf=influence(unim)
plot(inf)

##remove 8 effect size(Baggio et al. 2014PD-MCI) 

data=data[-8,]

#re-run multiple meta after remove influence study
mlm.inf <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm.inf


fmri.d=list(
  mlm=mlm,
  egger=egger,
  mlm.inf=mlm.inf
)

fmri[[metric]]=fmri.d
rm(list = setdiff(ls(), c("meta.d", "fmri","results", lsf.str())))

############## Assortativity######
metric="Assortativity"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm


#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim

inf=influence(unim)
plot(inf)

##remove 4 effect size(Sreenivasan et al. 2023 PD-FOG)

data=data[-3,]

#re-run multiple meta after remove influence study
mlm.inf <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm.inf


fmri.d=list(
  mlm=mlm,
  egger=egger,
  mlm.inf=mlm.inf
)

fmri[[metric]]=fmri.d
rm(list = setdiff(ls(), c("meta.d", "fmri","results", lsf.str())))

############## BC######
metric="BC"
data=meta.d[[metric]]

# meta
unim <- rma.uni(yi, vi, slab = Study, data = data, method = "REML")
unim


fmri.d=list(
  unim=unim
)

fmri[[metric]]=fmri.d
rm(list = setdiff(ls(), c("meta.d", "fmri","results", lsf.str())))
############## strength###### <-
metric="Strength"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm


#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim

inf=influence(unim)
plot(inf)


##Continuous variable
continuous_vars = c("PD.Age.years", "PD.Male.prop", "PD.Education.years", "OnsetAge", 
                    "Duration.years", "UPDRS_III", "H_Y", "MMSE","LEDD","Nodes")

continuous_results = cont_res(continuous_vars)


##Categorical variable
categorical_vars = c("Drug_code", "Scanning", "Template_code","Network_framework", "Threshold_code")
data=process_categorical_vars(categorical_vars, data)

categorical_results = cate_res(categorical_vars,data)


fmri.d=c(list(
  mlm=mlm,
  egger=egger
),continuous_results,categorical_results)

fmri[[metric]]=fmri.d
rm(list = setdiff(ls(), c("meta.d", "fmri","results", lsf.str())))

############## Degree######
metric="Degree"
data=meta.d[[metric]]

# meta
unim <- rma.uni(yi, vi, slab = Study, data = data, method = "REML")
unim


fmri.d=list(
  unim=unim
)

fmri[[metric]]=fmri.d
rm(list = setdiff(ls(), c("meta.d", "fmri","results", lsf.str())))
############## Participation coefficient######
metric="PC"
data=meta.d[[metric]]

# meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm

fmri.d=list(
  mlm=mlm
)

fmri[[metric]]=fmri.d
rm(list = setdiff(ls(), c("meta.d", "fmri","results", lsf.str())))

############## save data#########
results[["fMRI"]]=fmri
###############################################sMRI####################
modality="sMRI"
source("data_processing.R")
meta.d=meta_d(modality)
############## Clustering coefficient######
metric="Cp"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm

#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim

inf=influence(unim)
plot(inf)

##remove 5 effect size(Yan et al. 2024)
data=data[-5,]

#re-run multiple meta after remove influence study
mlm.inf <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm.inf

smri=list()
smri.d=list(
  mlm=mlm,
  egger=egger,
  mlm.inf=mlm.inf
)

smri[[metric]]=smri.d
rm(list = setdiff(ls(), c("meta.d", "smri","results", lsf.str())))

############## Local efficiency######
metric="Eloc"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm

#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim

inf=influence(unim)
plot(inf)

##remove 5 effect size(Yan  et al. 2024)
data=data[-5,]

#re-run multiple meta after remove influence study
mlm.inf <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm.inf

smri.d=list(
  mlm=mlm,
  egger=egger,
  mlm.inf=mlm.inf
)

smri[[metric]]=smri.d
rm(list = setdiff(ls(), c("meta.d", "smri","results", lsf.str())))

############## Characteristic path length######
metric="Lp"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm

#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim

inf=influence(unim)
plot(inf)

smri.d=list(
  mlm=mlm,
  egger=egger
)

smri[[metric]]=smri.d
rm(list = setdiff(ls(), c("meta.d", "smri","results", lsf.str())))
############## Global efficiency######
metric="Eglob"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm

#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim

inf=influence(unim)
plot(inf)


##remove 1 effect size(Suo et al. 2021a)
data=data[-1,]

#re-run multiple meta after remove influence study
mlm.inf <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm.inf

smri.d=list(
  mlm=mlm,
  egger=egger,
  mlm.inf=mlm.inf
)

smri[[metric]]=smri.d
rm(list = setdiff(ls(), c("meta.d", "smri","results", lsf.str())))
############## normalized clustering coefficient######
metric="Gamma"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm

#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim

inf=influence(unim)
plot(inf)

smri.d=list(
  mlm=mlm,
  egger=egger
)

smri[[metric]]=smri.d
rm(list = setdiff(ls(), c("meta.d", "smri","results", lsf.str())))
############## normalized Characteristic path length######
metric="Lambda"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm

#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim

inf=influence(unim)
plot(inf)

##remove 1 effect size(Suo et al. 2021a)

data=data[-1,]

#re-run multiple meta after remove influence study
mlm.inf <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm.inf

smri.d=list(
  mlm=mlm,
  egger=egger,
  mlm.inf=mlm.inf
)

smri[[metric]]=smri.d
rm(list = setdiff(ls(), c("meta.d", "smri","results", lsf.str())))
############## small-worldness######
metric="Sigma"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm

#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim

inf=influence(unim)
plot(inf)

smri.d=list(
  mlm=mlm,
  egger=egger
)

smri[[metric]]=smri.d
rm(list = setdiff(ls(), c("meta.d", "smri","results", lsf.str())))


############## save data#########
results[["sMRI"]]=smri
###############################################EEG####################
modality="EEG"
source("data_processing.R")
meta.d=meta_d(modality)
############## Clustering coefficient######
metric="Cp"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm


#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim

inf=influence(unim)
plot(inf)

##remove 7 effect size(Zhang et al. 2022)
data=data[-7,]

#re-run multiple meta after remove influence study
mlm.inf <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm.inf


###meta-regression
data=meta.d[[metric]]

##Continuous variable
continuous_vars = c("PD.Age.years", "PD.Male.prop", "PD.Education.years", "OnsetAge", 
                    "Duration.years", "UPDRS_III", "H_Y", "MMSE","LEDD","Nodes")

continuous_results = cont_res(continuous_vars)

categorical_vars = c("Drug", "Template_code", "Threshold_code")
data=process_categorical_vars(categorical_vars, data)

categorical_results = cate_res(categorical_vars,data)

EEG.d=c(list(
  mlm=mlm,
  egger=egger,
  mlm.inf=mlm.inf
),categorical_results)

EEG=list()
EEG[[metric]]=EEG.d
rm(list = setdiff(ls(), c("meta.d", "EEG","results", lsf.str())))
############## Local efficiency######
metric="Eloc"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm


#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim

inf=influence(unim)
plot(inf)

##remove 1 effect size(Peláez Suárez et al. 2021PD-NC)
data=data[-1,]

#re-run multiple meta after remove influence study
mlm.inf <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm.inf


EEG.d=list(
  mlm=mlm,
  egger=egger,
  mlm.inf=mlm.inf
)

EEG[[metric]]=EEG.d
rm(list = setdiff(ls(), c("meta.d", "EEG","results", lsf.str())))
############## Characteristic path length######
metric="Lp"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm


#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim

inf=influence(unim)
plot(inf)

categorical_vars = c("Drug", "Template_code", "Threshold_code")
data=process_categorical_vars(categorical_vars, data)

categorical_results = cate_res(categorical_vars,data)

EEG.d=c(list(
  mlm=mlm,
  egger=egger
),categorical_results)


EEG[[metric]]=EEG.d
rm(list = setdiff(ls(), c("meta.d", "EEG","results", lsf.str())))
############## Global efficiency######
metric="Eglob"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm

#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim

inf=influence(unim)
plot(inf)

##remove 4 effect size(Yin et al. 2023PD-NFOG)
data=data[-4,]

#re-run multiple meta after remove influence study
mlm.inf <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm.inf

EEG.d=list(
  mlm=mlm,
  egger=egger,
  mlm.inf=mlm.inf
)

EEG[[metric]]=EEG.d
rm(list = setdiff(ls(), c("meta.d", "EEG","results", lsf.str())))
############## Transitivity######
metric="Transitivity"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm

EEG.d=list(
  mlm=mlm
)

EEG[[metric]]=EEG.d
rm(list = setdiff(ls(), c("meta.d", "EEG","results", lsf.str())))
############## Modularity######
metric="Modularity"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm


EEG.d=list(
  mlm=mlm
)

EEG[[metric]]=EEG.d
rm(list = setdiff(ls(), c("meta.d", "EEG","results", lsf.str())))
############## Strength######
metric="Strength"
data=meta.d[[metric]]

# multiple meta
unim <- rma.uni(yi, vi, slab = Study, data = data)
unim

EEG.d=list(
  unim=unim
)

EEG[[metric]]=EEG.d
rm(list = setdiff(ls(), c("meta.d", "EEG","results", lsf.str())))
############## Degree######
metric="Degree"
data=meta.d[[metric]]

# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm


EEG.d=list(
  mlm=mlm
)

EEG[[metric]]=EEG.d
rm(list = setdiff(ls(), c("meta.d", "EEG","results", lsf.str())))
############## save data#########
results[["EEG"]]=EEG
##################single study of MEG/PET/rsfMRIWM####################
##########MEG
modality="MEG"
source("data_processing.R")
meta.d=meta_d(modality)
#Normalized clustering coefficient
metric="Gamma"
data=meta.d[[metric]]
unim <- rma.uni(yi, vi, slab = Study, data = data)
unim

MEG=list()
MEG.d=list(
  unim=unim
)

MEG[[metric]]=MEG.d
rm(list = setdiff(ls(), c("meta.d", "MEG","results", lsf.str())))
#Normalized characteristic path length
metric="Lambda"
data=meta.d[[metric]]
unim <- rma.uni(yi, vi, slab = Study, data = data)
unim

MEG.d=list(
  unim=unim
)

MEG[[metric]]=MEG.d
rm(list = setdiff(ls(), c("meta.d", "MEG","results", lsf.str())))
#save data
results[["MEG"]]=MEG


##########PET
modality="PET"
source("data_processing.R")
meta.d=meta_d(modality)
#Clustering coefficient
metric="Cp"
data=meta.d[[metric]]
unim <- rma.uni(yi, vi, slab = Study, data = data)
unim

PET.d=list(
  unim=unim
)
PET=list()
PET[[metric]]=PET.d
rm(list = setdiff(ls(), c("meta.d", "PET","results", lsf.str())))
#Characteristic path length
metric="Lp"
data=meta.d[[metric]]
unim <- rma.uni(yi, vi, slab = Study, data = data)
unim

PET.d=list(
  unim=unim
)

PET[[metric]]=PET.d
rm(list = setdiff(ls(), c("meta.d", "PET","results", lsf.str())))
#gamma
metric="Gamma"
data=meta.d[[metric]]
unim <- rma.uni(yi, vi, slab = Study, data = data)
unim

PET.d=list(
  unim=unim
)

PET[[metric]]=PET.d
rm(list = setdiff(ls(), c("meta.d", "PET","results", lsf.str())))
#lambda
metric="Lambda"
data=meta.d[[metric]]
unim <- rma.uni(yi, vi, slab = Study, data = data)
unim

PET.d=list(
  unim=unim
)

PET[[metric]]=PET.d
rm(list = setdiff(ls(), c("meta.d", "PET","results", lsf.str())))
#sigma
metric="Sigma"
data=meta.d[[metric]]
unim <- rma.uni(yi, vi, slab = Study, data = data)
unim

PET.d=list(
  unim=unim
)

PET[[metric]]=PET.d
rm(list = setdiff(ls(), c("meta.d", "PET","results", lsf.str())))
#eloc
metric="Eloc"
data=meta.d[[metric]]
unim <- rma.uni(yi, vi, slab = Study, data = data)
unim

PET.d=list(
  unim=unim
)

PET[[metric]]=PET.d
rm(list = setdiff(ls(), c("meta.d", "PET","results", lsf.str())))
#Edlob
metric="Eglob"
data=meta.d[[metric]]
unim <- rma.uni(yi, vi, slab = Study, data = data)
unim

PET.d=list(
  unim=unim
)

PET[[metric]]=PET.d
rm(list = setdiff(ls(), c("meta.d", "PET","results", lsf.str())))
#Modularity
metric="Modularity"
data=meta.d[[metric]]
unim <- rma.uni(yi, vi, slab = Study, data = data)
unim

PET.d=list(
  unim=unim
)

PET[[metric]]=PET.d
rm(list = setdiff(ls(), c("meta.d", "PET","results", lsf.str())))
#Assortativity
metric="Assortativity"
data=meta.d[[metric]]
unim <- rma.uni(yi, vi, slab = Study, data = data)
unim

PET.d=list(
  unim=unim
)

PET[[metric]]=PET.d
rm(list = setdiff(ls(), c("meta.d", "PET","results", lsf.str())))
#Hierarchy
metric="Hierarchy"
data=meta.d[[metric]]
unim <- rma.uni(yi, vi, slab = Study, data = data)
unim

PET.d=list(
  unim=unim
)

PET[[metric]]=PET.d
rm(list = setdiff(ls(), c("meta.d", "PET","results", lsf.str())))
#Synchronization
metric="Synchronization"
data=meta.d[[metric]]
unim <- rma.uni(yi, vi, slab = Study, data = data)
unim

PET.d=list(
  unim=unim
)

PET[[metric]]=PET.d
rm(list = setdiff(ls(), c("meta.d", "PET","results", lsf.str())))
#save data
results[["PET"]]=PET


##########rsfMRIWM
modality="rsfMRIWM"
source("data_processing.R")
meta.d=meta_d(modality)
#gamma
metric="Gamma"
data=meta.d[[metric]]
unim <- rma.uni(yi, vi, slab = Study, data = data)
unim

rsfMRIWM.d=list(
  unim=unim
)

rsfMRIWM=list()
rsfMRIWM[[metric]]=rsfMRIWM.d
rm(list = setdiff(ls(), c("meta.d", "rsfMRIWM","results", lsf.str())))
#lambda
metric="Lambda"
data=meta.d[[metric]]
unim <- rma.uni(yi, vi, slab = Study, data = data)
unim

rsfMRIWM.d=list(
  unim=unim
)

rsfMRIWM[[metric]]=rsfMRIWM.d
rm(list = setdiff(ls(), c("meta.d", "rsfMRIWM","results", lsf.str())))
#sigma
metric="Sigma"
data=meta.d[[metric]]
unim <- rma.uni(yi, vi, slab = Study, data = data)
unim

rsfMRIWM.d=list(
  unim=unim
)

rsfMRIWM[[metric]]=rsfMRIWM.d
rm(list = setdiff(ls(), c("meta.d", "rsfMRIWM","results", lsf.str())))
#eloc
metric="Eloc"
data=meta.d[[metric]]
unim <- rma.uni(yi, vi, slab = Study, data = data)
unim

rsfMRIWM.d=list(
  unim=unim
)

rsfMRIWM[[metric]]=rsfMRIWM.d
rm(list = setdiff(ls(), c("meta.d", "rsfMRIWM","results", lsf.str())))
#eloc
metric="Eglob"
data=meta.d[[metric]]
unim <- rma.uni(yi, vi, slab = Study, data = data)
unim

rsfMRIWM.d=list(
  unim=unim
)

rsfMRIWM[[metric]]=rsfMRIWM.d
rm(list = setdiff(ls(), c("meta.d", "rsfMRIWM","results", lsf.str())))
#Strength
metric="Strength"
data=meta.d[[metric]]
unim <- rma.uni(yi, vi, slab = Study, data = data)
unim

rsfMRIWM.d=list(
  unim=unim
)

rsfMRIWM[[metric]]=rsfMRIWM.d
rm(list = setdiff(ls(), c("meta.d", "rsfMRIWM","results", lsf.str())))
#save data
results[["rsfMRIWM"]]=rsfMRIWM

###############################################Cog-fMRI#####################
modality="fMRI"
source("data_processing.R")
meta.d=meta_d(modality)

#######cognitive normal#####
metric="Cp"
data=meta.d[[metric]]

data=data %>% filter(Non_motor=="Cog",Symptom_code=="withoutsymptom")
# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm
#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger
#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim
inf=influence(unim)
plot(inf)
##remove 3 effect size(Lopes et al. 2017Group 1)
data=data[-3,]
#re-run multiple meta after remove influence study
mlm.inf <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm.inf

Cog.d=list(mlm=mlm,
           egger=egger,
           mlm.inf=mlm.inf)

Cog=list()
Cog[[metric]]=Cog.d
rm(list = setdiff(ls(), c("meta.d", "Cog","results", lsf.str())))


metric="Eloc"
data=meta.d[[metric]]

data=data %>% filter(Non_motor=="Cog",Symptom_code=="withoutsymptom")
# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm
#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger
#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim
inf=influence(unim)
plot(inf)


Cog.d=list(mlm=mlm,
           egger=egger)

Cog[[metric]]=Cog.d
rm(list = setdiff(ls(), c("meta.d", "Cog","results", lsf.str())))

metric="Lp"
data=meta.d[[metric]]

data=data %>% filter(Non_motor=="Cog",Symptom_code=="withoutsymptom")
# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm
#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger
#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim
inf=influence(unim)
plot(inf)

Cog.d=list(mlm=mlm,
           egger=egger)

Cog[[metric]]=Cog.d
rm(list = setdiff(ls(), c("meta.d", "Cog","results", lsf.str())))

metric="Eglob"
data=meta.d[[metric]]

data=data %>% filter(Non_motor=="Cog",Symptom_code=="withoutsymptom")
# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm
#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger
#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim
inf=influence(unim)
plot(inf)

Cog.d=list(mlm=mlm,
           egger=egger)

Cog[[metric]]=Cog.d
rm(list = setdiff(ls(), c("meta.d", "Cog","results", lsf.str())))
results[["Cognitivenormal"]]=Cog

#######Cognitive impairment#####
metric="Cp"
data=meta.d[[metric]]

data=data %>% filter(Non_motor=="Cog",Symptom_code=="withsymptom")
# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm
#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger
#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim
inf=influence(unim)
plot(inf)
##remove 4 effect size(Lopes et al. 2017Group 4)
data=data[-4,]
#re-run multiple meta after remove influence study
mlm.inf <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm.inf

Cog.d=list(mlm=mlm,
           egger=egger,
           mlm.inf=mlm.inf)

Cog=list()
Cog[[metric]]=Cog.d
rm(list = setdiff(ls(), c("meta.d", "Cog","results", lsf.str())))


metric="Eloc"
data=meta.d[[metric]]

data=data %>% filter(Non_motor=="Cog",Symptom_code=="withsymptom")
# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm
#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger
#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim
inf=influence(unim)
plot(inf)
##remove 5 effect size(Suo et al. 2022PD-MCI)
data=data[-5,]
#re-run multiple meta after remove influence study
mlm.inf <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm.inf

Cog.d=list(mlm=mlm,
           egger=egger,
           mlm.inf=mlm.inf)

Cog[[metric]]=Cog.d
rm(list = setdiff(ls(), c("meta.d", "Cog","results", lsf.str())))

metric="Lp"
data=meta.d[[metric]]

data=data %>% filter(Non_motor=="Cog",Symptom_code=="withsymptom")
# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm
#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger
#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim
inf=influence(unim)
plot(inf)
##remove 8 effect size(Suo et al. 2022PD-MCI)
data=data[-8,]
#re-run multiple meta after remove influence study
mlm.inf <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm.inf

Cog.d=list(mlm=mlm,
           egger=egger,
           mlm.inf=mlm.inf)

Cog[[metric]]=Cog.d
rm(list = setdiff(ls(), c("meta.d", "Cog","results", lsf.str())))

metric="Eglob"
data=meta.d[[metric]]

data=data %>% filter(Non_motor=="Cog",Symptom_code=="withsymptom")
# multiple meta
mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm
#publication bias by egger meta-regression
egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger
#influence test
unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim
inf=influence(unim)
plot(inf)

Cog.d=list(mlm=mlm,
           egger=egger)

Cog[[metric]]=Cog.d
rm(list = setdiff(ls(), c("meta.d", "Cog","results", lsf.str())))
results[["Cognitiveimpairment"]]=Cog

#######Cognitive impairment vs cognitive normal #####
modality="fMRI"
meta.d=meta_d(modality,"symptom")

metric="Cp"
data=meta.d[[metric]]
data=data %>% filter(data$Non_motor=="Cog")

mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm

egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim
inf=influence(unim)
plot(inf)

Cog.d=list(mlm=mlm,
           egger=egger)

Cog=list()
Cog[[metric]]=Cog.d
rm(list = setdiff(ls(), c("meta.d", "Cog","results", lsf.str())))

metric="Eloc"
data=meta.d[[metric]]
data=data %>% filter(data$Non_motor=="Cog")

mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm

egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim
inf=influence(unim)
plot(inf)

##remove 8 effect size(Lopes et al. 2017Group 2)
data=data[-2,]
#re-run multiple meta after remove influence study
mlm.inf <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm.inf

Cog.d=list(mlm=mlm,
           egger=egger,
           mlm.inf=mlm.inf)

Cog[[metric]]=Cog.d
rm(list = setdiff(ls(), c("meta.d", "Cog","results", lsf.str())))

metric="Lp"
data=meta.d[[metric]]
data=data %>% filter(data$Non_motor=="Cog")

mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm

egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim
inf=influence(unim)
plot(inf)
##remove 6 effect size(	Lopes et al. 2017Group 4)
data=data[-6,]
#re-run multiple meta after remove influence study
mlm.inf <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm.inf

Cog.d=list(mlm=mlm,
           egger=egger,
           mlm.inf=mlm.inf)

Cog[[metric]]=Cog.d
rm(list = setdiff(ls(), c("meta.d", "Cog","results", lsf.str())))

metric="Eglob"
data=meta.d[[metric]]
data=data %>% filter(data$Non_motor=="Cog")

mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm

egger <- rma.mv(yi, vi, mod =  ~ sqrt(vi), random = ~ 1|study.id/es.id, data = data, test="t",method = "REML")
egger

unim <- rma.uni(yi, vi, slab = Study, data = data, test="t", method = "REML")
unim
inf=influence(unim)
plot(inf)

Cog.d=list(mlm=mlm,
           egger=egger)

Cog[[metric]]=Cog.d
rm(list = setdiff(ls(), c("meta.d", "Cog","results", lsf.str())))
results[["Cognitiveimpirmentvsnormal"]]=Cog


############## robust_process###########

process_robust <- function(model) {
  if (inherits(model, "rma.mv")) {
    if (model$s.nlevels[1] > 1) {
      return(robust(model, cluster = study.id, clubSandwich = TRUE)) 
    } else  return(model)  
  } else  return(model) }

results.r <- list()

results.r <- map(results, ~map(.x, ~map(.x, process_robust)))


#########save results robust data############
save(results, results.r, file = "results.r.RData")



