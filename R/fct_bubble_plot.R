
#' bubble_plot 
#'
#' @description Plot a ggplot2 bubble plot for the selected Key Events
#'
#' @return Returns ggplot2 object.
#'
#' @import dplyr
#' @import ggplot2
#' @import scales
#' @export
#' @noRd
bubble_plot <- function(optimal_models_stats, Enrichment_data, deg_statistics, df_tissue_specific, rel_ke, PFAS, plot_title=paste("PFAS", PFAS, sep=" - ")) {
  # gene logFC - max abs for each gene logfc across PFAS, time point 
  deg_top_by_group <- deg_statistics %>%
    dplyr::group_by(Experiment, time, Feature) %>%
    dplyr::slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Exp=paste(Experiment, time, sep = "_"))

  df_tissue_specific = merge(AOPfingerprintR::Biological_system_annotations, df_tissue_specific, by.x="key_event_name", by.y="KeyEvent")
  df_tissue_specific = unique(df_tissue_specific[,c("key_event_name","level")])

  optimal_models_stats$Experiment = paste(optimal_models_stats$Experiment,optimal_models_stats$time, sep = "_")

  ED <- Enrichment_data

  #MW_PFOS = 500.126	#PFOS molecular mass (g/mol)
  #MW_PFOA = 414.07	#PFOA molecular mass (g/mol)
  MW_vec <- c("PFOS"=500.126, "PFOA"=414.07)
  MW <- MW_vec[PFAS] |> unname()

  #idx = intersect(which(ED$key_event_name %in% rel_ke), grep(pattern = "PFOS",x = ED$Experiment))
  #liver_paths = unique(ED[idx,"Ke_description"])
  #df = unique(ED[ED$Ke_description %in% liver_paths, c("Experiment","Genes","BMD","BMDL","BMDU","Ke_description")])

  DF <- ED |> dplyr::filter(key_event_name %in% rel_ke & stringr::str_detect(Experiment, pattern=PFAS)) |> dplyr::select(Experiment, Genes, BMD, BMDL, BMDU, Ke_description) 

  #df$Day = as.numeric(unlist(lapply(strsplit(df$Experiment,split = "_"), function(elem)elem[2])))
  DF$Day = DF$Experiment |> stringr::str_replace(pattern=".*_", replacement="") |> as.numeric()
  DF$BMD = as.numeric(DF$BMD) * MW
  DF$BMDL = as.numeric(DF$BMDL) * MW
  DF$BMDU = as.numeric(DF$BMDU) * MW

  DF = merge(DF, df_tissue_specific, by.x = "Ke_description", by.y = "key_event_name")
  #DF <- DF %>% mutate(Ke_description = fct_relevel(Ke_description, unique(Ke_description[order(level)])))
  DF <- DF %>% mutate(Ke_description=Ke_description |> as.vector(), Ke_description = factor(Ke_description, unique(Ke_description[order(level)])))

  dd = c()
  for(i in 1:nrow(DF)){
    di = optimal_models_stats[optimal_models_stats$Exp == DF$Experiment[i],]
    gi = unlist(strsplit(DF$Genes[i],split = ";"))
    dd_i = di[di$Feature %in% gi,c("Feature","BMD","Experiment","Adverse_direction")]
    deg_i = deg_top_by_group[deg_top_by_group$Exp == DF$Experiment[i], ]
    deg_i = deg_i[deg_i$Feature %in% gi,c("Feature","log2FoldChange")]

    dd_i = merge(dd_i, deg_i, by.x = "Feature",by.y = "Feature", all.x = T)
    dd = rbind(dd, dd_i)
  }

  dd$BMD = dd$BMD * MW

  #dd$Time = unlist(lapply(strsplit(dd$Experiment,split = "_"), function(elem)elem[2]))
  #dd$Time = factor(dd$Time, levels = c(1,4, 10, 14))
  dd$Time = dd$Experiment |> stringr::str_replace(pattern=".*_", replacement="")
  time_levels <- unique(dd$Time) |> as.numeric() |> sort() |> as.character()
  dd$Time=factor(dd$Time, levels = time_levels)

  dd <- dd %>% dplyr::mutate(BMD_inv = max(BMD, na.rm = TRUE) - BMD)

  p5 <- ggplot(dd, aes(x = Time, y = Feature)) +
    geom_point(aes(size = BMD_inv, color = log2FoldChange), alpha = 0.9) +
    scale_color_gradient2(
      low = "darkblue", mid = "white", high = "firebrick",
      midpoint = 0, name = "log2FC"
    ) +
    # Area-scaled bubbles
    scale_size_area(
      max_size = 8,
      name = "BMD (µg/L)",
      breaks = pretty_breaks(4),
      labels = function(br) {
        # show legend labels in terms of the original BMD
        inv <- max(dd$BMD, na.rm = TRUE) - br
        # optional: round nicely
        format(inv, trim = TRUE, digits = 3)
      }
    ) +
    theme_minimal(base_size = 16) +
    labs(
      title = plot_title,
      x = "In vitro Time (Day)",
      y = "Genes"
    ) + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 16),
      plot.title = element_text(size = 16, face = "bold")
    )
    return(p5)
}

