

process_categorical_vars <- function(categorical_vars, data) {
  
  # Convert all categorical variables to factors
  data[categorical_vars] <- lapply(data[categorical_vars], as.factor)
  
  # Report the subgroup table for each variable
  lapply(categorical_vars, function(var) {
    cat("Subgroup for", var, ":\n")
    print(table(data[[var]]))
  })
  
  # Report variables with more than 1 subgroup
  lapply(categorical_vars, function(var) {
    num_groups <- length(unique(data[[var]]))
    if (num_groups > 1) {
      cat("More than 1 subgroup for", var, "\n")
    }
  })
  # Return the modified data
  return(data)
}

cate_res <- function(categorical_vars, data) {
  
  # Initialize the final results list
  final_results <- list()
  
  lapply(categorical_vars, function(var) {
    
    
    data[[var]] <- factor(data[[var]])
    subgroups <- levels(data[[var]])
    subgroup_counts <- table(data[[var]]) 
    
    if (length(subgroups) < 2 || any(subgroup_counts < 3)) {
      message(paste("Skipping", var, ": Insufficient subgroups or subgroup size <3"))
      return(NULL)
    }
    
    # Main effect analysis
    main_effect <- tryCatch({
      rma.mv(yi, vi, slab = Study,
             random = ~ 1 | study.id/es.id,
             mods = as.formula(paste("~", var)),
             data = data, test = "t", method = "REML")
    }, error = function(e) {
      message(paste("Main effect optimizer failed for", var, "- switching to Nelder-Mead"))
      rma.mv(yi, vi, slab = Study,
             random = ~ 1 | study.id/es.id,
             mods = as.formula(paste("~", var)),
             data = data, test = "t", method = "REML",
             control = list(optimizer = "Nelder-Mead"))
    })
    
    # Store the main effect in final results
    main_name <- paste0(var)
    final_results[[main_name]] <<- main_effect
    
    # Subgroup analysis
    lapply(subgroups, function(level) {
      subset_data <- data %>% 
        dplyr::filter(.data[[var]] == !!level) 
      n_effects <- nrow(subset_data)
      
      # Subgroup model creation
      subgroup_model <- if (n_effects == 1) {
        tryCatch(
          rma.uni(yi, vi, slab = Study, data = subset_data),
          error = function(e) NULL
        )
      } else if (n_effects >= 2) {
        tryCatch({
          rma.mv(yi, vi, slab = Study,
                 random = ~ 1 | study.id/es.id,
                 data = subset_data, test = "t", method = "REML")
        }, error = function(e) {
          rma.mv(yi, vi, slab = Study,
                 random = ~ 1 | study.id/es.id,
                 data = subset_data, test = "t", method = "REML",
                 control = list(optimizer = "Nelder-Mead")) 
        })
      } else NULL
      
      # Store valid subgroup in final results
      if (!is.null(subgroup_model)) {
        final_results[[as.character(level)]] <<- subgroup_model
      }
    })
  })
  
  return(final_results)
}

cont_res <- function(continuous_vars) {
  continuous_results <- lapply(continuous_vars, function(var) {
    
    # Check the number of valid effect sizes for the continuous variable
    valid_study_count = sum(!is.na(data[[var]]))
    
    if (valid_study_count >= 10) {
      result <- tryCatch({
        # Perform meta-analysis for the continuous variable
        rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, 
               mods = as.formula(paste("~ data$", var)), data = data, test = "t", method = "REML")
      }, error = function(e) {
        # If error occurs, use a different optimizer and report it
        message(paste("Optimizer (nlminb) did not achieve convergence for", var, "- switching to Nelder-Mead optimizer"))
        rma.mv(yi, vi, slab = Study, random = ~ 1 | study.id/es.id, 
               mods = as.formula(paste("~ data$", var)), data = data, test = "t", 
               method = "REML", control = list(optimizer = "Nelder-Mead"))
      })
      return(result)
    } else {
      return(NULL)
    }
  })
  
  # Name the results list by the continuous variable names
  names(continuous_results) = continuous_vars
  
  # Filter out NULL results (variables that don't meet the condition)
  continuous_results = Filter(Negate(is.null), continuous_results)
  
  return(continuous_results)
}

regp <- function(metric, modality,cont){
  
  model=results.r[[modality]][[metric]][[cont]]
  data=model$data%>% mutate(n = PD.N + HC.N,ci.lb=yi-1.96*sqrt(vi),
                            ci.ub=yi+1.96*sqrt(vi)) %>%   filter(!is.na(.data[[cont]])) #%>% rename(mod=all_of(cont))
  reg_data <- data %>%select(Study,all_of(cont),yi,ci.lb,ci.ub,n)
  
  YLAB <- paste0(metric_names[metric])
  XLAB <- paste0(cont)
  
  cont_min=min(data[[cont]], na.rm = TRUE)
  cont_max=max(data[[cont]], na.rm = TRUE)
  range=cont_max-cont_min
  
  mod=paste0("data$",cont)
  
  at_list <- list()
  at_list[[mod]] <- seq(cont_min - range/10, cont_max + range/10, range/200)
  
  #emm <- emmprep(model, at=list("data$LEDD"=seq(cont_min - range/10, cont_max + range/10, range/200)))
  model$vb <- as.matrix(model$vb)
  emm <- emmprep(model, at=at_list)
  emms <- emmeans(emm, specs=cont, weights="equal") %>% as.data.frame() %>%
    rename(yi=emmean, ci.lb=lower.CL, ci.ub=upper.CL)
  
  category_colors <- c(
    "dMRI" =      "#d92b03", 
    "fMRI" =      "#214ea7", 
    "sMRI" =      "#7552a7", 
    "rsfMRIWM" =  "#f38b2f", 
    "PET" =       "#f1b534", 
    "EEG" =       "#088c00", 
    "MEG" =       "#a5405e"
  )
  
  in.nbreaks=c(10,50,250)
  textsize=11
  
  p <- ggplot(data=emms, aes(x=!!sym(cont), y=yi))+
    #geom_ribbon(aes(ymin=ci.lb, ymax=ci.ub), fill="gray90", color="transparent")+
    geom_ribbon(aes(ymin=yi-SE, ymax=yi+SE), fill="gray90", color="transparent")+
    geom_line(linewidth=1)+
    geom_point(data=data, aes(size=log(n)), shape=21, color="white", fill=category_colors[modality], alpha=0.7)+
    scale_size_identity(breaks=log(in.nbreaks), labels=in.nbreaks, limits=log(c(1,500)),
                        guide=guide_legend(override.aes=list(shape=16, color="gray25")))+
    labs( size="N", x=XLAB, y=YLAB)+
    theme(panel.background=element_blank(),
          panel.grid=element_blank(),
          plot.tag=element_text(size=textsize*0.8, face="bold", color="black"),
          text=element_text(size=textsize, color="black"),
          axis.line.x=element_line(lineend="square"),
          axis.ticks.x=element_line(lineend="square"),
          axis.text.x=element_text(color="black"),
          axis.title.x=element_text(size=textsize*0.8, color="black"),
          axis.line.y=element_line(lineend="square"),
          axis.ticks.y=element_line(lineend="square"),
          # axis.title.y = element_blank(),
          axis.title.y=element_text(size=textsize*0.8, color="black"),
          legend.position="none",
          legend.key=element_blank(),
          legend.text=element_text(color="black"),
          #legend.title=element_text(size=textsize*0.8, face="bold"),
          legend.margin=margin(0,0,0,0),
          plot.margin=margin(3,3,3,3))+ labs(title = paste0(metric_names[metric],"-",modality))
  return(p)
}

# sensitivity change data extracted, e.g. guanweighted data to binary data, shangbanatomical atlas to functional parcellation
sens_process_metric <- function(metric, params, study) {
  data <- meta.d[[metric]]
  
  # Update the data with the given study name and params
  data[data$Study == study, c("pd.m", "pd.sd", "hc.m", "hc.sd")] <- 
    list(params$pd.m, params$pd.sd, params$hc.m, params$hc.sd)
  
  # Recalculate effect size
  data[data$Study == study, ] <- escalc(
    measure = "SMD", 
    n1i = PD.N, n2i = HC.N,
    m1i = pd.m, m2i = hc.m,
    sd1i = pd.sd, sd2i = hc.sd,
    data = data[data$Study == study, ]
  )
  
  # Run meta-analysis model
  mlm <- rma.mv(
    yi, vi,
    random = ~ 1 | study.id/es.id,
    data = data,
    test = "t", 
    method = "REML",
    slab = Study
  )
  
  # Store results
  sensresult[[sens]][[modality]][[metric]] <- mlm
  return(sensresult)
}


#####heterogeneity
heterogeneity <- function(mlm) {
  if (inherits(mlm, "rma")) {
    W <- diag(1 / mlm$vi)  
    X <- model.matrix(mlm)  
    P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W  
    return(data.frame(
      Tau = sqrt(sum(mlm$sigma2)),
      I2 = 100 * sum(mlm$sigma2) / (sum(mlm$sigma2) + (mlm$k - mlm$p) / sum(diag(P))),  
      Q = mlm$QE,  
      Qdf = mlm$QEdf,  
      Qp = mlm$QEp
    ))
  }}

