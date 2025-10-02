metrics <- c("Cp", "Eloc","Lp", "Eglob", "Gamma", "Lambda", "Sigma", "Modularity")

####removing task-study##### 
sens="removingtask"
modality="fMRI"
source("data_processing.R")
meta.d=meta_d(modality)

metric="Cp"

data=meta.d[[metric]]
data$Modality
data=data %>% filter(Modality=="rs-fMRI")


mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm
sensresult=list()

sensresult[[sens]][[modality]][[metric]]=mlm

metric="Strength"

data=meta.d[[metric]]
data$Modality
data=data %>% filter(Modality=="rs-fMRI")


mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm
sensresult[[sens]][[modality]][[metric]]=mlm


modality="EEG"
meta.d=meta_d(modality)

metric="Cp"

data=meta.d[[metric]]
data$Modality
data=data %>% filter(Modality=="EEG")


mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm
sensresult[[sens]][[modality]][[metric]]=mlm

metric="Lp"

data=meta.d[[metric]]
data$Modality
data=data %>% filter(Modality=="EEG")

mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm
sensresult[[sens]][[modality]][[metric]]=mlm

metric="Eglob"

data=meta.d[[metric]]
data$Modality
data=data %>% filter(Modality=="EEG")

mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
mlm
sensresult[[sens]][[modality]][[metric]]=mlm

##########removingmultiplesparsity##########
sens="removingmultiplesparsity"
modality="fMRI"
source("data_processing.R")
meta.d=meta_d(modality)

metrics <- c("Cp", "Eloc","Lp", "Eglob", "Gamma", "Lambda", "Sigma", "Modularity")

for (metric in metrics) {
  data=meta.d[[metric]]
  data$pd.m
  data=data %>% filter(!is.na(pd.m))
  mlm <- rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, data = data, test = "t", method = "REML")
  mlm
  sensresult[[sens]][[modality]][[metric]]=mlm
}






####guan binary data####
sens="guanbinary"
####dmri
modality="dMRI"
meta.d=meta_d(modality)
# 
metrics <- list(
  Cp = list(pd.m = 0.091, pd.sd = 0.003, hc.m = 0.091, hc.sd = 0.003),
  Lp = list(pd.m = 0.38, pd.sd = 0.008, hc.m = 0.378, hc.sd = 0.005),
  Eloc = list(pd.m = 0.137, pd.sd = 0.003, hc.m = 0.136, hc.sd = 0.003),
  Eglob = list(pd.m = 0.106, pd.sd = 0.002, hc.m = 0.106, hc.sd = 0.001)
)

for (metric in names(metrics)) {
  sensresult <- sens_process_metric(metric, metrics[[metric]], study = "Guan et al. 2019")
}

#######fmri
modality="fMRI"
meta.d=meta_d(modality)
##
# Define metrics to process and corresponding parameters
metrics <- list(
  Cp = list(pd.m = 0.278, pd.sd = 0.01, hc.m = 0.286, hc.sd = 0.011),
  Lp = list(pd.m = 0.891, pd.sd = 0.050, hc.m = 0.928, hc.sd = 0.066),
  Eloc = list(pd.m = 0.347, pd.sd = 0.008, hc.m = 0.349, hc.sd = 0.010),
  Eglob = list(pd.m = 0.254, pd.sd = 0.007, hc.m = 0.250, hc.sd = 0.009)
)


# Batch process all metrics (with study parameter)
for (metric in names(metrics)) {
  sensresult <- sens_process_metric(metric, metrics[[metric]], study = "Guan et al. 2019")
}
####shangb functional atlas data####
sens="shangbfunction"
####dmri
modality="dMRI"
meta.d=meta_d(modality)
# 
metrics <- list(
  Cp = list(pd.m = 0.02, pd.sd = 0.01, hc.m = 0.03, hc.sd = 0.01),
  Lp = list(pd.m = 2.65, pd.sd = 0.91, hc.m = 1.91, hc.sd = 0.44),
  Eloc = list(pd.m = 0.02, pd.sd = 0.01, hc.m = 0.03, hc.sd = 0.01),
  Eglob = list(pd.m = 0.02, pd.sd = 0.01, hc.m = 0.03, hc.sd = 0.01),
  Gamma = list(pd.m = 1.96, pd.sd = 0.29, hc.m = 1.77, hc.sd = 0.34),
  Lambda = list(pd.m = 0.32, pd.sd = 0.02, hc.m = 0.31, hc.sd = 0.01),
  Sigma = list(pd.m = 1.38, pd.sd = 0.26, hc.m = 1.56, hc.sd = 0.23)
)

for (metric in names(metrics)) {
  sensresult <- sens_process_metric(metric, metrics[[metric]], study = "Shang et al. 2023b")
}

#######fmri
modality="fMRI"
meta.d=meta_d(modality)
##
# Define metrics to process and corresponding parameters
metrics <- list(
  Cp =     list(pd.m = 0.13, pd.sd = 0.01, hc.m = 0.14, hc.sd = 0.01),
  Lp =     list(pd.m = 0.48, pd.sd = 0.02, hc.m = 0.50, hc.sd = 0.03),
  Eloc =   list(pd.m = 0.19, pd.sd = 0.01, hc.m = 0.19, hc.sd = 0.01),
  Eglob =  list(pd.m = 0.14, pd.sd = 0.01, hc.m = 0.13, hc.sd = 0.01),
  Gamma =  list(pd.m = 0.54, pd.sd = 0.08, hc.m = 0.48, hc.sd = 0.10),
  Lambda = list(pd.m = 0.27, pd.sd = 0.01, hc.m = 0.28, hc.sd = 0.01),
  Sigma =  list(pd.m = 0.49, pd.sd = 0.07, hc.m = 0.43, hc.sd = 0.10)
)


# Batch process all metrics (with study parameter)
for (metric in names(metrics)) {
  sensresult <- sens_process_metric(metric, metrics[[metric]], study = "Shang et al. 2023b")
}



#####robust#####
process_robust <- function(model) {
  return(robust(model, cluster = study.id, clubSandwich = TRUE)) 
}

sensresult.r <- list()

sensresult.r <- map(sensresult, ~map(.x, ~map(.x, process_robust)))

#
save(sensresult, sensresult.r, file = "sensresult.r.RData")
