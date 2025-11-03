meta_d <- function(modality, symptom = NULL) {
  # -------------------------------
  # Function: Perform data preprocessing and effect size calculation
  # Inputs:
  #   modality: one of "fMRI", "dMRI", "sMRI", "EEG", "MEG", "rsfMRIWM"(rs-fMRI white matter), "PET"
  #   symptom: default NULL; if not NULL, perform symptom subgroup analysis
  # Output:
  #   A list of data frames, one per metric, containing effect sizes
  # -------------------------------
  
  # 1. File path mapping
  # File path management
  data_paths <- list(
    fMRI = "F:\\meta\\PD-GTA-Meta\\Data extract\\PD\\inputdata_fMRI.xlsx",
    dMRI = "F:\\meta\\PD-GTA-Meta\\Data extract\\PD\\inputdata_dMRI.xlsx",
    sMRI = "F:\\meta\\PD-GTA-Meta\\Data extract\\PD\\inputdata_sMRI.xlsx",
    EEG = "F:\\meta\\PD-GTA-Meta\\Data extract\\PD\\inputdata_EEG.xlsx",
    MEG = "F:\\meta\\PD-GTA-Meta\\Data extract\\PD\\inputdata_MEG.xlsx",
    rsfMRIWM = "F:\\meta\\PD-GTA-Meta\\Data extract\\PD\\inputdata_rsfMRIWM.xlsx",
    PET = "F:\\meta\\PD-GTA-Meta\\Data extract\\PD\\inputdata_PET.xlsx"
  )
  
  # Select file path based on modality
  data.path <- data_paths[[modality]]
  d <- read_excel(path = data.path)
  

  # Define metrics and labels (combined into one vector)
  if (modality == "fMRI") {
    metrics <- c("Cp", "Lp", "Gamma", "Lambda", "Sigma", "Eloc", "Eglob", "Transitivity", "Modularity", "Assortativity","Strength","Degree","PC")
  } else if (modality == "dMRI"){
    metrics <- c("Cp", "Lp", "Gamma", "Lambda", "Sigma", "Eloc", "Eglob","Transitivity", "Modularity","Assortativity","BC", "Strength","Degree","Density","Eccentricity")
  } else if (modality == "sMRI"){
    metrics <- c("Cp", "Lp", "Gamma", "Lambda", "Sigma", "Eloc", "Eglob")
  } else if (modality == "EEG"){
    metrics <- c("Cp", "Lp", "Sigma", "Eloc", "Eglob","Transitivity","Modularity", "Strength","Degree")
  } else if (modality == "MEG"){
    metrics <- c("Gamma", "Lambda" )
  } else if (modality == "rsfMRIWM"){
    metrics <- c("Gamma", "Lambda", "Sigma", "Eloc", "Eglob", "Strength")
  } else if (modality == "PET"){
    metrics <- c("Cp", "Lp", "Gamma", "Lambda", "Sigma", "Eloc", "Eglob", "Modularity", "Assortativity","Hierarchy","Synchronization")
  }

 
  # Function to process each metric
  process_metric <- function(data, metric) {
    # Filter data, mutate, select relevant columns and rename them
    data_processed <- data %>%
      filter(!is.na(!!sym(paste0(metric, ".PD")))) %>%
      mutate(GTA = metric) %>%
      select(Study:N, starts_with(metric), GTA) %>%
      rename_with(~ c("pd.m", "pd.sd", "hc.m", "hc.sd"), ends_with(c("PD", "PD.SD", "HC", "HC.SD")))
    
    # Calculate effect size using escalc and assign it to the data
    data_processed <- escalc(
      measure = "SMD",
      n1i = PD.N, n2i = HC.N,
      m1i = pd.m, m2i = hc.m,
      sd1i = pd.sd, sd2i = hc.sd,
      data = data_processed
    )
    
    return(data_processed) 
  }
  
  
  if(!is.null(symptom)){
  
  metrics=c("Cp", "Lp", "Eloc", "Eglob","Gamma","Lambda", "Sigma")
  
  # Function to process each metric
  process_metric <- function(data, metric) {
    # Filter data, mutate, select relevant columns and rename them
    data_processed <- data %>%
      filter(!is.na(Symptom)) %>%
      filter(!is.na(!!sym(paste0(metric, ".PD")))) %>%
      mutate(GTA = metric) %>%
      select(Study:N, starts_with(metric), GTA) %>%
      rename_with(~ c("pd.m", "pd.sd", "hc.m", "hc.sd"), ends_with(c("PD", "PD.SD", "HC", "HC.SD")))
  
    data_processed <- data_processed  %>%
      group_by(Study = factor(Study, levels = unique(Study))) %>%
      mutate(study.id = cur_group_id())
  
    data_processed <- data_processed %>%
      group_by(study.id) %>%
      mutate(
        HC.N = if (all(Symptom_code == "withsymptom") || all(Symptom_code == "withoutsymptom")) HC.N else ifelse(Symptom_code == "withsymptom", PD.N[Symptom_code == "withoutsymptom"], HC.N),
        hc.m = if (all(Symptom_code == "withsymptom") || all(Symptom_code == "withoutsymptom")) hc.m else ifelse(Symptom_code == "withsymptom", pd.m[Symptom_code == "withoutsymptom"], hc.m),
        hc.sd = if (all(Symptom_code == "withsymptom") || all(Symptom_code == "withoutsymptom")) hc.sd else ifelse(Symptom_code == "withsymptom", pd.sd[Symptom_code == "withoutsymptom"], hc.sd)
      ) %>%filter(Symptom_code != "withoutsymptom") %>%
      ungroup()
  
    # Calculate effect size using escalc and assign it to the data
    data_processed <- escalc(
      measure = "SMD",
      n1i = PD.N, n2i = HC.N,
      m1i = pd.m, m2i = hc.m,
      sd1i = pd.sd, sd2i = hc.sd,
      data = data_processed
    )
    return(data_processed)
  }
  }

  
  
  # Process all metrics using lapply
  metric_data <- lapply(metrics, function(metric) process_metric(d, metric))
  
  # Assign names to the processed data
  names(metric_data) <- metrics
  
  
  if (modality == "fMRI" & is.null(symptom)) {
    
    dm <- read_excel(path = data.path, sheet = "multi-sparsity")
    
    fMRI_metrics <- c("Cp", "Lp", "Gamma", "Lambda", "Sigma", "Eloc", "Eglob", "Modularity", "BC")
    
    for (metric in fMRI_metrics) {
      metric_yi <- paste0(metric, ".yi")
      metric_vi <- paste0(metric, ".vi")

        dm_metric <- dm %>%
          filter(!is.na(!!sym(metric_yi))) %>%  
          mutate(pd.m = N, pd.sd = N, hc.m = N, hc.sd = N) %>%  
          mutate(GTA = metric) %>%  
          select(Study:N, pd.m:hc.sd, GTA, !!sym(metric_yi):!!sym(metric_vi)) %>%
          rename(yi = !!sym(metric_yi), vi = !!sym(metric_vi))
        
        if (metric %in% metrics) {
          metric_data[[metric]] <- bind_rows(metric_data[[metric]], dm_metric)
        } else {
          metric_data[[metric]] <- dm_metric
        }
      }
    }
  
  
  if (modality == "fMRI" & !is.null(symptom)) {
    
    dm <- read_excel(path = data.path, sheet = "multi-sparsity-symptom")
    
    fMRI_metrics <- c("Cp", "Lp", "Eglob")
    
    for (metric in fMRI_metrics) {
      metric_yi <- paste0(metric, ".yi")
      metric_vi <- paste0(metric, ".vi")
      
      dm_metric <- dm %>%
        filter(!is.na(!!sym(metric_yi))) %>%  
        mutate(pd.m = N, pd.sd = N, hc.m = N, hc.sd = N) %>%  
        mutate(GTA = metric) %>%  
        select(Study:N, pd.m:hc.sd, GTA, !!sym(metric_yi):!!sym(metric_vi)) %>%
        rename(yi = !!sym(metric_yi), vi = !!sym(metric_vi))
      
      if (metric %in% metrics) {
        metric_data[[metric]] <- bind_rows(metric_data[[metric]], dm_metric)
      } else {
        metric_data[[metric]] <- dm_metric
      }
    }
  }
  
  # Function to add group IDs for study and sample
  add_group_ids <- function(data) {
    data %>%
      group_by(Study = factor(Study, levels = unique(Study))) %>%
      mutate(study.id = cur_group_id()) %>%
      ungroup() %>%
      mutate(es.id = row_number())
  }
  

  metric_data <- lapply(metric_data, add_group_ids)
  
  return(metric_data)
}



