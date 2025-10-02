

e.cont <- function(ed) {
  
  if (nrow(ed$b) == 2) {
    k <- calculate_k(ed$s.nlevels)
    out <- data.frame(
       k=k,
      B = round(ed$b[2],3),
      `95% CI` = paste0(round(ed$ci.lb[2], 3), " to ", round(ed$ci.ub[2], 3)),
      se = round(ed$se[2],3),
      t = round(ed$zval[2],3),
      Pt = round(ed$pval[2],3),
      check.names = FALSE
    )
  }  else {
    stop("Error: Unsupported number of rows in ed$b")
  }
  return(out)
}

e.cate <- function(ed) {
  # Check if ed$b has enough rows
  n_rows <- nrow(ed$b)
  
  if (n_rows < 2) {
    stop("Error: ed$b must have at least 2 rows")
  }
  
  # Determine the rows to process (from 2nd row to the last row)
  rows_to_process <- 2:n_rows
  
  # Use lapply to process all required rows
  out <- do.call(rbind, lapply(rows_to_process, function(i) {
    k <- calculate_k(ed$s.nlevels)
    data.frame(
      k = k,
      B = round(ed$b[i], 3),
      `95% CI` = paste0(round(ed$ci.lb[i], 3), " to ", round(ed$ci.ub[i], 3)),
      se = round(ed$se[i], 3),
      t = round(ed$zval[i], 3),
      Pt = round(ed$pval[i], 3),
      check.names = FALSE
    )
  }))
  
  # Set row names
  row.names(out) <- rownames(ed$b)[rows_to_process]
  
  return(out)
}

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
calculate_k <- function(s.nlevels) {
  
  if (length(s.nlevels) == 3) {  
    k <- paste0(s.nlevels[1], "-", s.nlevels[2], "-", s.nlevels[3])
  } else if (length(s.nlevels) == 2) {
    k <- paste0(s.nlevels[1], "-", s.nlevels[2])
  } else {
    stop("Unexpected length of s.nlevels. It should be 2 or 3.")
  }
  return(k)
}

e.mlm <- function(ed) {
  k <- calculate_k(ed$s.nlevels)
  if (nrow(ed$b) == 1) {
    ed.h <- heterogeneity(ed)
    out <- data.frame(
      k = k,
      g = round(ed$b[1], 3),
      `95% CI` = paste0(round(ed$ci.lb, 3), " to ", round(ed$ci.ub, 3)),
      se= round(ed$se, 3),
      t = round(ed$zval, 3),
      Pt =round(ed$pval, 3),
      Q = round(ed$QE, 3),
      PQ = round(ed$QEp, 3),
      `I² %` = round(ed.h$I2, 3),
      Tau = round(ed.h$Tau, 3),
      check.names = FALSE
    )
  }  else {
    stop("Error: Unsupported number of rows in ed$b")
  }
  
  return(out)
}

e.unim <- function(ed) {
  if (nrow(ed$b) == 1) {
    out <- data.frame(
      k = c("1-1"),
      g = round(ed$b[1], 3),
      `95% CI` = paste0(round(ed$ci.lb, 3), " to ", round(ed$ci.ub, 3)),
      se= round(ed$se, 3),
      t = round(ed$zval, 3),
      Pt = round(ed$pval, 3),
      check.names = FALSE
    )
  }  else {
    stop("Error: Unsupported number of rows in ed$b")
  }
  
  return(out)
}
e.main<- function(ed){
  if("s.nlevels" %in% names(ed)){
    out <- e.mlm(ed)
  }else{out <- e.unim(ed)}
  return(out)
}

e.egger <- function(ed) {
  k <- calculate_k(ed$s.nlevels)
  if (nrow(ed$b) == 2) {
    out <- data.frame(
      k = k,
      intercept = round(ed$b[1], 3),
      B = round(ed$b[2], 3),
      se = round(ed$se[2], 3),
      t = round(ed$zval[2], 3),
      Pt = round(ed$pval[2], 3)
    )
  }  else {
    stop("Error: Unsupported number of rows in ed$b")
  }
  return(out)
}




format3_df <- function(df) {
  cols <- names(df)[(which(names(df)=="k")+1):ncol(df)]
  df %>% mutate(across(all_of(cols), ~ 
                         sapply(as.character(.), function(x) {
                           if (grepl("to",x)) paste(sprintf("%.3f",as.numeric(strsplit(x,"\\s+to\\s+")[[1]])),collapse=" to ")
                           else if (grepl("^-?\\d",x)) sprintf("%.3f",as.numeric(x)) 
                           else x
                         })))}



calculate_stats <- function(df) {
  k <- paste0(length(unique(df$Study)), "-", nrow(df))
  
  # PD.N
  pd_n_total <- sum(df$PD.N)
  
  # Age
  age_mean <- round(mean(df$PD.Age, na.rm = TRUE), 2)
  age_sd <- round(sd(df$PD.Age, na.rm = TRUE), 2)
  age_na <- sum(is.na(df$PD.Age))
  
  # Gender
  gender_mean <- round(mean(df$PD.Male.prop, na.rm = TRUE), 2)
  gender_sd <- round(sd(df$PD.Male.prop, na.rm = TRUE), 2)
  gender_na <- sum(is.na(df$PD.Male.prop))
  
  # Education
  Education_mean <- round(mean(df$Education, na.rm = TRUE), 2)
  Education_sd <- round(sd(df$Education, na.rm = TRUE), 2)
  Education_na <- sum(is.na(df$Education))
  
  # Duration
  Duration_mean <- round(mean(df$Duration, na.rm = TRUE), 2)
  Duration_sd <- round(sd(df$Duration, na.rm = TRUE), 2)
  Duration_na <- sum(is.na(df$Duration))
  
  # UPDRS
  UPDRS_mean <- round(mean(df$`UPDRS-III`, na.rm = TRUE), 2)
  UPDRS_sd <- round(sd(df$`UPDRS-III`, na.rm = TRUE), 2)
  UPDRS_na <- sum(is.na(df$`UPDRS-III`))
  
  # H&Y
  HY_mean <- round(mean(df$`H&Y`, na.rm = TRUE), 2)
  HY_sd <- round(sd(df$`H&Y`, na.rm = TRUE), 2)
  HY_na <- sum(is.na(df$`H&Y`))
  
  # MMSE
  MMSE_mean <- round(mean(df$MMSE, na.rm = TRUE), 2)
  MMSE_sd <- round(sd(df$MMSE, na.rm = TRUE), 2)
  MMSE_na <- sum(is.na(df$MMSE))
  
  # Drug
  on.state <- table(df$Drug)[["On-state"]]
  off.state <- table(df$Drug)[["Drug-naïve/off-state"]]
  drug_na <- sum(is.na(df$Drug))
  
  # LEDD
  LEDD_mean <- round(mean(df$LEDD, na.rm = TRUE), 2)
  LEDD_sd <- round(sd(df$LEDD, na.rm = TRUE), 2)
  LEDD_na <- sum(is.na(df$LEDD))
  
  
  data.frame(
    Characteristic = c(
      "No. of studies/samples",
      "No. of PD",
      "Age", "",
      "Gender", "",
      "Education", "",
      "Duration", "",
      "UPDRS-III", "",
      "H&Y", "",
      "MMSE", "",
      "Drug", "", "",
      "LEDD", ""
    ),
    Metrics = c("Numbers", "Numbers",
                "Mean (SD)", "Missing samples",
                "Male% Mean (SD) ", "Missing samples",
                "Mean (SD)", "Missing samples",
                "Mean (SD)", "Missing samples",
                "Mean (SD)", "Missing samples",
                "Mean (SD)", "Missing samples",
                "On-state", "Drug-naïve/ off-state", "Missing samples",
                "Mean (SD)", "Missing samples"),
    modality = c(
      k, pd_n_total,
      sprintf("%.2f (%.2f)", age_mean, age_sd),
      age_na,
      sprintf("%.2f (%.2f)", gender_mean, gender_sd),
      gender_na,
      sprintf("%.2f (%.2f)", Education_mean, Education_sd),
      Education_na,
      sprintf("%.2f (%.2f)", Duration_mean, Duration_sd),
      Duration_na,
      sprintf("%.2f (%.2f)", UPDRS_mean, UPDRS_sd),
      UPDRS_na,
      sprintf("%.2f (%.2f)", HY_mean, HY_sd),
      HY_na,
      sprintf("%.2f (%.2f)", MMSE_mean, MMSE_sd),
      MMSE_na,
      on.state, off.state, drug_na,
      sprintf("%.2f (%.2f)", LEDD_mean, LEDD_sd),
      LEDD_na
    )
  )
}
