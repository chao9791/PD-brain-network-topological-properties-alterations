library(patchwork)
library(emmeans)

load("results.r.RData")


modalities <- c("dMRI", "fMRI", "sMRI","EEG")
singlestudy <- c("rsfMRIWM","PET","MEG")

metrics <- c("Cp", "Eloc","Lp", "Eglob", "Gamma", "Lambda", "Sigma", "Modularity")
metric_names <- c(
  "Cp" = "Clustering Coefficient",
  "Eloc" = "Local Efficiency",
  "Lp" = "Characteristic Path Length", 
  "Eglob" = "Global Efficiency",
  "Gamma" = "Normalized Clustering Coefficient",
  "Lambda" = "Normalized Characteristic Path Length",
  "Sigma" = "Small-Worldness",
  "Modularity" = "Modularity",
  "Strength" = "Strength"
)

##### Report modalities and graph metrics with points exceeding coordinate limits
orchard <- list()
for (metric in metrics) {
  
  d <- list()
  for (modality in modalities) {
    if(!modality %in% singlestudy){
      d[[modality]] <- results.r[[modality]][[metric]]$mlm
    } else {
      d[[modality]] <- results.r[[modality]][[metric]]$unim
    }
  }
  
  pooles <- map_df(names(d), ~ data.frame(
    modalities = .x,
    yi = d[[.x]]$b,
    se = d[[.x]]$se
  ))
  
  es <- list()
  for (modality in names(d)) {
    es[[modality]] <- d[[modality]]$data
  }
  
  es.all <- bind_rows(
    lapply(names(es), function(modality) {
      es[[modality]] %>%
        mutate(n = PD.N + HC.N) %>%       # Calculate n
        select(Study, n, yi) %>%           # Select needed columns
        mutate(modalities = modality)      # Add modality column
    })
  )
  
  # Set factor levels to reversed unique modalities levels
  common_levels <- rev(unique(c(pooles$modalities, es.all$modalities)))
  pooles$modalities <- factor(pooles$modalities, levels = common_levels)
  es.all$modalities <- factor(es.all$modalities, levels = common_levels)
  
  # Combine both data frames
  pooles$source <- "pooles"
  es.all$source <- "es.all"
  orchard.data <- bind_rows(pooles, es.all)
  
  category_colors <- c(
    "dMRI" =      "#d92b03", 
    "fMRI" =      "#214ea7", 
    "sMRI" =      "#7552a7", 
    "rsfMRIWM" =  "#f38b2f", 
    "PET" =       "#f1b534", 
    "EEG" =       "#088c00", 
    "MEG" =       "#a5405e"
  )
  
  # Identify out-of-bound points before plotting
  out_points <- es.all %>% 
    filter(yi < -2.5 | yi > 2.5) %>% 
    group_by(modalities) %>%
    summarise(
      n_out = n(),
      studies = toString(unique(Study))
    )
  
  if(nrow(out_points) > 0) {
    message("\n==========================================")
    message(paste("Metric:", metric, "has points outside [-2.5,2.5] range in:"))
    for(i in 1:nrow(out_points)) {
      message(paste0("  - ", out_points$modalities[i], 
                     " (", out_points$n_out[i], " points from studies: ", 
                     out_points$studies[i], ")"))
    }
    message("==========================================")
  }
  
  # Create plot (original code unchanged)
  orchard[[metric]] <- ggplot(data = orchard.data, aes(y = modalities, x = yi)) +
    geom_beeswarm(data = es.all, aes(x = yi, y = modalities, size = log(n), fill = modalities), 
                  shape = 21, color = "white", alpha = 0.6, cex = 2.5,
                  corral = "wrap", corral.width = 0.5) +
    geom_vline(xintercept = 0, lineend = "square", linetype = 3) +
    geom_errorbarh(data = pooles, aes(xmin = yi-se, xmax = yi+se), linewidth = 0.75, lineend = "square", height = 0) +
    geom_point(data = pooles, aes(x = yi, fill = modalities), size = 2, shape = 21, stroke = 1, color = "black") +
    scale_size_identity(breaks = log(c(10,50,250)), 
                        labels = c(10,50,250), 
                        limits = log(c(1,500)),
                        guide = guide_legend(override.aes = list(shape = 16, color = "gray50", alpha = 1))) +
    scale_fill_manual(values = category_colors[modalities]) +
    scale_x_continuous(limits = c(-2.5, 2.5), breaks = c(-2, -1, 0, 1, 2)) +
    labs(size = "N", color = "modalities", x = "Hedges' g") +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.tag = element_text(size = 8, face = "bold", color = "black"),  # modified here
      text = element_text(size = 9, color = "black"),  # modified here
      axis.line.x = element_line(lineend = "square"),
      axis.text.x = element_text(color = "black"),  # modified here
      axis.title.x = element_text(size = 8, color = "black"),  # modified here
      axis.line.y = element_line(lineend = "square"),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(hjust = 1, color = "black"),  # modified here
      axis.title.y = element_blank(),
      # legend.position = "none",
      legend.key = element_blank(),
      legend.text = element_text(color = "black"),  # modified here
      legend.title = element_text(size = 8, face = "bold", color = "black"),  # modified here
      legend.margin = margin(0,0,0,0),
      plot.margin = margin(3,3,3,3) # added title font settings
    ) + labs(title = metric_names[metric])
  print(orchard[[metric]])
}


#### Combine 6 plots into one
wrap_plots(orchard[[1]],orchard[[2]],orchard[[3]],orchard[[4]],orchard[[5]],orchard[[6]],orchard[[7]],orchard[[8]], ncol = 2)


######forest#####
if (!dir.exists("forest")) {dir.create("forest")}

for (modality in modalities) {
  for (metric in metrics) {
    model=results.r[[modality]][[metric]][[1]]
    if (inherits(model, "rma.mv")) {
      k=model$s.nlevels[2]
      if ( k>5) {
        file_path <- file.path("forest", paste0(metric, "-", modality,"-forest", ".png"))
        png(file_path , width=23, height <- max(12, 0.7 * k), units="cm", res=300, bg="white")
        p=forest(model, ilab=cbind(model$data$Subgroup),shade = "zebra", order = "obs") 
        dev.off()
      }
    }
  }
}


###### funnel#####
if (!dir.exists("funnel")) {dir.create("funnel")}

for (modality in modalities) {
  for (metric in metrics) {
    model=results.r[[modality]][[metric]][[1]]
    if (inherits(model, "rma.mv")) {
      k=model$s.nlevels[2]
      if ( k>5) {
        file_path <- file.path("funnel", paste0(metric, "-", modality,"-funnel", ".png"))
        png(file_path , width=15, height=12, units="cm", res=300, bg="white")
        p=funnel(model)
        dev.off()
        
      }
    }
  }
}



######regression-plot####

if (!dir.exists("regression")) {dir.create("regression")}

regplist=list()

contvars=list(
  contvar1=c("Eloc","dMRI","LEDD"),
  contvar2=c("Sigma","dMRI","Nodes"),
  contvar3=c("Eglob","fMRI","PD.Education.years"),
  contvar4=c("Strength","fMRI","PD.Male.prop"),
  contvar5=c("Lambda","fMRI","UPDRS_III")
)


for (contvar in contvars) {
  print(contvar)
  file_path <- file.path("regression", paste0(contvar[1], "-", contvar[2],"-",contvar[3],"-regression", ".png"))
  png(file_path , width=15, height=12, units="cm", res=300, bg="white")
  regplist[[contvar[1]]] =regp(metric=contvar[1],modality=contvar[2],cont=contvar[3])
  print(regplist[[contvar[1]]])
  dev.off()
}


catevars=list(
  catevar1=c("Cp","fMRI","1.5T","3T"),
  catevar2=c("Cp","fMRI","Macroanatomy","Function"),
  catevar3=c("Eloc","fMRI","nonProportionalthreshold","Proportionalthreshold"),
  catevar4=c("Lambda","fMRI","drug.naïve.or.off.state","on.state"))

orchard <- list()
for (catevar in catevars) {
  
  cates <- tail(catevar, 2)  
  
  d <- list()
  for (cate in cates) {
    d[[cate]]=results.r[[catevar[2]]][[catevar[1]]][[cate]]
  }
  
  pooles <- map_df(names(d), ~ data.frame(
    cates = .x,
    yi = d[[.x]]$b,
    se = d[[.x]]$se
  ))
  
  es <- list()
  for (cate in names(d)) {
    es[[cate]] <- d[[cate]]$data
  }
  
  es.all <- bind_rows(
    lapply(names(es), function(cate) {
      es[[cate]] %>%
        mutate(n = PD.N + HC.N) %>%       # Calculate n
        select(Study, n, yi) %>%           # Select needed columns
        mutate(cates = cate)      # Add modality column
    })
  )
  
  # Combine both data frames
  pooles$source <- "pooles"
  es.all$source <- "es.all"
  orchard.data <- bind_rows(pooles, es.all)
  
  
  if (all(c("Function", "Macroanatomy") %in% orchard.data$cates)) {
    # Set custom order: Function on top, Macroanatomy on bottom
    orchard.data$cates <- factor(
      orchard.data$cates,
      levels = c("Function", "Macroanatomy")
    )
  } else {
    # Keep original order for other groups
    orchard.data$cates <- factor(
      orchard.data$cates,
      levels = sort(unique(orchard.data$cates))
    )
  }
  
  category_colors <- c(
    "1.5T"               = "#88c0c0",  
    "3T"                 = "#0e8585", 
    "Macroanatomy"       = "#d4a854",  
    "Function"           = "#7e4909",  
    "drug.naïve.or.off.state" = "#a3cb87",  
    "on.state"           = "#1e803d",  
    "nonProportionalthreshold" = "#c19ac1",  
    "Proportionalthreshold"    = "#830783"   
  )
  
  # Identify out-of-bound points before plotting
  out_points <- es.all %>% 
    filter(yi < -2.5 | yi > 2.5) %>% 
    group_by(cates) %>%
    summarise(
      n_out = n(),
      studies = toString(unique(Study))
    )
  
  if(nrow(out_points) > 0) {
    message("\n==========================================")
    message(paste("Metric:", catevar[1], "has points outside [-2.5,2.5] range in:"))
    for(i in 1:nrow(out_points)) {
      message(paste0("  - ", out_points$cates[i], 
                     " (", out_points$n_out[i], " points from studies: ", 
                     out_points$studies[i], ")"))
    }
    message("==========================================")
  }
  
  if (catevar[2] == "fMRI" && "Function" %in% orchard.data$cates) {
    orchard.data$cates <- factor(orchard.data$cates,
                                 levels = c( "Macroanatomy","Function"))
  }
  # Create plot (original code unchanged)
  p <- ggplot(data = orchard.data, aes(y = cates, x = yi)) +
    geom_beeswarm(data = es.all, aes(x = yi, y = cates, size = log(n), fill = cates), 
                  shape = 21, color = "white", alpha = 0.6, cex = 2.5,
                  corral = "wrap", corral.width = 0.5) +
    geom_vline(xintercept = 0, lineend = "square", linetype = 3) +
    geom_errorbarh(data = pooles, aes(xmin = yi-se, xmax = yi+se), linewidth = 0.75, lineend = "square", height = 0) +
    geom_point(data = pooles, aes(x = yi, fill = cates), size = 2, shape = 21, stroke = 1, color = "black") +
    scale_size_identity(breaks = log(c(10,50,250)), 
                        labels = c(10,50,250), 
                        limits = log(c(1,500)),
                        guide = "none") +
    scale_fill_manual(values = category_colors[cates],
                      guide = guide_legend(
                        title = "cates",
                        title.position = "top",
                        reverse = TRUE  
                      )) +
    scale_x_continuous(limits = c(-2.5, 2.5), breaks = c(-2, -1, 0, 1, 2)) +
    labs(size = "N", color = "cates", x = "Hedges' g") +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.tag = element_text(size = 8, face = "bold", color = "black"),  # modified here
      text = element_text(size = 9, color = "black"),  # modified here
      axis.line.x = element_line(lineend = "square"),
      axis.text.x = element_text(color = "black"),  # modified here
      axis.title.x = element_text(size = 8, color = "black"),  # modified here
      axis.line.y = element_line(lineend = "square"),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      # legend settings (transparent effect)
      legend.background = element_rect(fill = NA, color = NA),  
      legend.key = element_rect(fill = NA, color = NA),        
      legend.title = element_blank(),
      legend.position = "inside",
      legend.position.inside = if (catevar[1] %in% c("Cp","Eloc")) {
        c(0.02, 0.2) } else {  c(0.02, 0.9)   },
      legend.justification = c(0, if (catevar[1] %in% c("Cp","Eloc")) 0 else 1),
      legend.box.margin = margin(0, 0, 0, 0),  
      legend.text = element_text(color = "black"),  # modified here
      legend.margin = margin(0,0,0,0),
      plot.margin = margin(3,3,3,3) # added title font settings
    ) + labs(title = paste0(metric_names[catevar[1]],"-",catevar[2]), face = "bold")
  print(orchard[[catevar[4]]])
  
  # Only enforce y-axis order if both Macroanatomy and Function are present
  if (all(c("Macroanatomy", "Function") %in% unique(orchard.data$cates))) {
    p <- p + scale_y_discrete(limits = c("Macroanatomy", "Function"))
  }
  
  orchard[[catevar[4]]] <- p
  print(orchard[[catevar[4]]])
}



wrap_plots(regplist[[1]],regplist[[2]],regplist[[3]],regplist[[4]],regplist[[5]],orchard[[1]],orchard[[2]],orchard[[3]],orchard[[4]], ncol = 3)


print(regplist[[3]])

dev.new()  # open new graphic device
print(regplist[[3]])
orchard[[1]]

