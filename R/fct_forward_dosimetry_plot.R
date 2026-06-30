
#' Plot PBK forward-dosimetry simulations across scenarios and chemicals
#'
#' Builds a PBK plot from one or more forward-dosimetry runs (e.g., the
#' outputs you store as `out2`, `out32`). Each run can contain one or two
#' chemicals (typically `PFOS` and/or `PFOA`) under a given exposure scenario.
#' The function combines runs, computes uncertainty bands when needed, and
#' produces a faceted plot (one facet per chemical) showing median
#' concentrations with 5th–95th percentile error bars across time.
#'
#' @param ... One or more forward-dosimetry run objects. Each run is expected
#'   to be a list containing per-chemical entries (e.g., `out$PFOS`, `out$PFOA`)
#'   with a component `bma$liver_df`. The `liver_df` can be either:
#'   \itemize{
#'     \item a data frame with columns \code{p50}, \code{p5}, \code{p95}, or
#'     \item a matrix/data frame of posterior samples by time where the function
#'           will compute the 5th/50th/95th percentiles by row.
#'   }
#' @param scenario_labels Character vector of scenario labels, one per run
#'   passed in \code{...}. If \code{NULL}, the function will try to use the
#'   names of the runs or fall back to \code{"Scenario 1"}, \code{"Scenario 2"}, etc.
#' @param chemicals Character vector of chemical names to extract from each run.
#'   Defaults to \code{c("PFOS", "PFOA")}. Only chemicals present in a given run
#'   are included.
#' @param y_label Character string for the y-axis label. Default is
#'   \code{"Concentrations Liver (µg/L)"}.
#'
#' @details
#' The function is designed to handle the common cases:
#' \itemize{
#'   \item 1 chemical × 1 scenario
#'   \item 1 chemical × multiple scenarios
#'   \item 2 chemicals × 1 scenario
#'   \item 2 chemicals × multiple scenarios
#' }
#' For each run and chemical, it looks for \code{out[[chem]]$bma$liver_df}.
#' If \code{liver_df} already contains \code{p50}, \code{p5}, \code{p95}, these
#' are used directly. Otherwise, it assumes rows correspond to time points and
#' columns to posterior samples, and computes the 5th/50th/95th percentiles
#' across samples for each time row.
#'
#' The resulting plot uses points, lines, and error bars, with a facet per
#' chemical (free y-scales). Scenario is mapped to color.
#'
#' @return A \link[ggplot2]{ggplot} object.
#'
#' @section Assumptions:
#' \itemize{
#'   \item Each run corresponds to a single exposure scenario.
#'   \item Time is indexed as \code{Year = 1:n} where \code{n} is the number of
#'         rows in \code{liver_df} (or the length of \code{p50} when present).
#'   \item When \code{liver_df} is not in \code{p50/p5/p95} format, rows are
#'         treated as time points and quantiles are computed across columns.
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1: One chemical, two scenarios
#' p1 <- plot_PBK_simulations(
#'   out2, out32,
#'   scenario_labels = c("130 ng/kg/day", "32 ng/kg/day"),
#'   chemicals = "PFOS"
#' )
#' p1
#'
#' # Example 2: Two chemicals in one scenario (assuming both exist in out2)
#' p2 <- plot_PBK_simulations(
#'   out2,
#'   scenario_labels = "130 ng/kg/day",
#'   chemicals = c("PFOS", "PFOA")
#' )
#' p2
#'
#' # Example 3: Two chemicals, two scenarios
#' p3 <- plot_PBK_simulations(
#'   out2, out32,
#'   scenario_labels = c("130 ng/kg/day", "32 ng/kg/day"),
#'   chemicals = c("PFOS", "PFOA")
#' )
#' p3
#' }
#'
#' @import ggplot2
#' @import dplyr
#' @export
plot_PBK_simulations <- function(...,
                          scenario_labels = NULL,
                          chemicals = c("PFOS", "PFOA"),
                          y_label = "Concentrations Liver (µg/L)",
                          return_data = FALSE,
                          scale_values = "fixed") {
  runs <- list(...)
  
  # If labels not provided, try names of ... or default
  if (is.null(scenario_labels)) {
    scenario_labels <- names(runs)
    if (is.null(scenario_labels) || any(scenario_labels == "")) {
      scenario_labels <- paste0("Scenario ", seq_along(runs))
    }
  }
  
  stopifnot(length(runs) == length(scenario_labels))
  
  all_dfs <- list()
  
  for (i in seq_along(runs)) {
    out <- runs[[i]]
    scen <- scenario_labels[[i]]
    
    for (chem in chemicals) {
      # Only proceed if this chemical exists in this run
      if (!is.null(out[[chem]]) && !is.null(out[[chem]]$bma$liver_df)) {
        liver <- out[[chem]]$bma$liver_df
        
        # Two lightweight cases:
        # 1) Already has p50/p5/p95 columns
        # 2) Is a matrix/data.frame of samples → compute 5/50/95% quantiles by row
        if (is.data.frame(liver) && all(c("p50", "p5", "p95") %in% names(liver))) {
          df <- data.frame(
            Year           = seq_len(length(liver$p50)),
            Concentration  = liver$p50,
            ConcentrationL = liver$p5,
            ConcentrationU = liver$p95,
            Scenario       = scen,
            Chemical       = chem,
            stringsAsFactors = FALSE
          )
        } else {
          # Assume 'liver' is samples-by-time or time-by-samples (your example was time-by-samples)
          #q <- t(apply(liver, 1, function(row) quantile(row, probs = c(0.05, 0.5, 0.95))))
          q <- t(apply(liver |> dplyr::select(-time), 1, function(row) quantile(row, probs = c(0.05, 0.5, 0.95))))
          q <- as.data.frame(q)
          colnames(q) <- c("ConcentrationL", "Concentration", "ConcentrationU")
          df <- cbind(
            Year     = seq_len(nrow(q)),
            q,
            Scenario = scen,
            Chemical = chem
          )
          df <- as.data.frame(df)
        }
        
        all_dfs[[length(all_dfs) + 1]] <- df
      }
    }
  }
  
  df_PBK <- dplyr::bind_rows(all_dfs) %>%
    dplyr::mutate(
      Scenario = factor(Scenario, levels = scenario_labels),
      Chemical = factor(Chemical, levels = chemicals)
    )

  if(return_data){
    return(df_PBK)
  }
  
  # Title: if one chemical, use it; if two, join with " & "
  # chem_title <- paste(unique(as.character(df_PBK$Chemical)), collapse = " & ")
  
  PBK_plot <- ggplot(df_PBK, aes(x = Year, y = Concentration, color = Scenario, group = Scenario)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = ConcentrationL, ymax = ConcentrationU), width = 0.3, color = "gray40") +
    geom_line(size = 1) +
    facet_wrap(~ Chemical, scales = scale_values) +
    theme_minimal(base_size = 16) +
    labs(
      # title = chem_title,
      x = "Year",
      y = y_label,
      color = "Scenario"
    ) +
    theme(
      plot.title     = element_text(size = 16, face = "bold"),
      axis.title     = element_text(size = 16),
      axis.text      = element_text(size = 16),
      legend.position = "bottom",
      legend.title   = element_text(size = 14),
      legend.text    = element_text(size = 12)
    )
  
  return(PBK_plot)
}

#' Builds a PBK plot from one or more forward-dosimetry runs (e.g., the
#' outputs you store as `out2`, `out32`). Each run can contain one or two
#' chemicals (typically `PFOS` and/or `PFOA`) under a given exposure scenario.
#'
#' @param ... One or more forward-dosimetry run objects. Each run is expected
#'   to be a list containing per-chemical entries (e.g., `out$PFOS`, `out$PFOA`)
#'   with a component `bma$liver_df`. The `liver_df` can be either:
#'   \itemize{
#'     \item a data frame with columns \code{p50}, \code{p5}, \code{p95}, or
#'     \item a matrix/data frame of posterior samples by time where the function
#'           will compute the 5th/50th/95th percentiles by row.
#'   }
#' @param scenario_labels Character vector of scenario labels, one per run
#'   passed in \code{...}. If \code{NULL}, the function will try to use the
#'   names of the runs or fall back to \code{"Scenario 1"}, \code{"Scenario 2"}, etc.
#' @param chemicals Character vector of chemical names to extract from each run.
#'   Defaults to \code{c("PFOS", "PFOA")}. Only chemicals present in a given run
#'   are included.
#' @param y_label Character string for the y-axis label. Default is
#'   \code{"Concentrations Liver (µg/L)"}.
#'
#' @details
#' The function is designed to handle the common cases:
#' \itemize{
#'   \item 1 chemical × 1 scenario
#'   \item 1 chemical × multiple scenarios
#'   \item 2 chemicals × 1 scenario
#'   \item 2 chemicals × multiple scenarios
#' }
#' For each run and chemical, it looks for \code{out[[chem]]$bma$liver_df}.
#' If \code{liver_df} already contains \code{p50}, \code{p5}, \code{p95}, these
#' are used directly. Otherwise, it assumes rows correspond to time points and
#' columns to posterior samples, and computes the 5th/50th/95th percentiles
#' across samples for each time row.
#'
#' The resulting plot uses points, lines, and error bars, with a facet per
#' chemical (free y-scales). Scenario is mapped to color.
#'
#' @return A \link[ggplot2]{ggplot} object.
#'
#' @section Assumptions:
#' \itemize{
#'   \item Each run corresponds to a single exposure scenario.
#'   \item Time is indexed as \code{Year = 1:n} where \code{n} is the number of
#'         rows in \code{liver_df} (or the length of \code{p50} when present).
#'   \item When \code{liver_df} is not in \code{p50/p5/p95} format, rows are
#'         treated as time points and quantiles are computed across columns.
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1: One chemical, two scenarios
#' p1 <- plot_PBK_simulations_serum(
#'   out2, out32,
#'   scenario_labels = c("130 ng/kg/day", "32 ng/kg/day"),
#'   chemicals = "PFOS"
#' )
#' p1
#'
#' # Example 2: Two chemicals in one scenario (assuming both exist in out2)
#' p2 <- plot_PBK_simulations_serum(
#'   out2,
#'   scenario_labels = "130 ng/kg/day",
#'   chemicals = c("PFOS", "PFOA")
#' )
#' p2
#'
#' # Example 3: Two chemicals, two scenarios
#' p3 <- plot_PBK_simulations_serum(
#'   out2, out32,
#'   scenario_labels = c("130 ng/kg/day", "32 ng/kg/day"),
#'   chemicals = c("PFOS", "PFOA")
#' )
#' p3
#' }
#'
#' @import ggplot2
#' @import dplyr
#' @export
plot_PBK_simulations_serum <- function(...,
                                 scenario_labels = NULL,
                                 chemicals = c("PFOS", "PFOA"),
                                 y_label = "Concentrations Serum (µg/L)",
                                 return_data = FALSE,
                                 scale_values = "fixed") {
  runs <- list(...)
  
  # If labels not provided, try names of ... or default
  if (is.null(scenario_labels)) {
    scenario_labels <- names(runs)
    if (is.null(scenario_labels) || any(scenario_labels == "")) {
      scenario_labels <- paste0("Scenario ", seq_along(runs))
    }
  }
  
  stopifnot(length(runs) == length(scenario_labels))
  
  all_dfs <- list()
  
  for (i in seq_along(runs)) {
    out <- runs[[i]]
    scen <- scenario_labels[[i]]
    
    for (chem in chemicals) {
      # Only proceed if this chemical exists in this run
      if (!is.null(out[[chem]]) && !is.null(out[[chem]]$bma$serum_df)) {
        serum <- out[[chem]]$bma$serum_df
        
        # Two lightweight cases:
        # 1) Already has p50/p5/p95 columns
        # 2) Is a matrix/data.frame of samples → compute 5/50/95% quantiles by row
        if (is.data.frame(serum) && all(c("p50", "p5", "p95") %in% names(serum))) {
          df <- data.frame(
            Year           = seq_len(length(serum$p50)),
            Concentration  = serum$p50,
            ConcentrationL = serum$p5,
            ConcentrationU = serum$p95,
            Scenario       = scen,
            Chemical       = chem,
            stringsAsFactors = FALSE
          )
        } else {
          # Assume 'liver' is samples-by-time or time-by-samples (your example was time-by-samples)
          #q <- t(apply(serum, 1, function(row) quantile(row, probs = c(0.05, 0.5, 0.95))))
          q <- t(apply(serum |> dplyr::select(-time), 1, function(row) quantile(row, probs = c(0.05, 0.5, 0.95))))
          q <- as.data.frame(q)
          colnames(q) <- c("ConcentrationL", "Concentration", "ConcentrationU")
          df <- cbind(
            Year     = seq_len(nrow(q)),
            q,
            Scenario = scen,
            Chemical = chem
          )
          df <- as.data.frame(df)
        }
        
        all_dfs[[length(all_dfs) + 1]] <- df
      }
    }
  }
  
  df_PBK <- dplyr::bind_rows(all_dfs) %>%
    dplyr::mutate(
      Scenario = factor(Scenario, levels = scenario_labels),
      Chemical = factor(Chemical, levels = chemicals)
    )
  
  if(return_data){
    return(df_PBK)
  }

  # Title: if one chemical, use it; if two, join with " & "
  # chem_title <- paste(unique(as.character(df_PBK$Chemical)), collapse = " & ")
  
  PBK_plot <- ggplot(df_PBK, aes(x = Year, y = Concentration, color = Scenario, group = Scenario)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = ConcentrationL, ymax = ConcentrationU), width = 0.3, color = "gray40") +
    geom_line(size = 1) +
    facet_wrap(~ Chemical, scales = scale_values) +
    theme_minimal(base_size = 16) +
    labs(
      # title = chem_title,
      x = "Year",
      y = y_label,
      color = "Scenario"
    ) +
    theme(
      plot.title     = element_text(size = 16, face = "bold"),
      axis.title     = element_text(size = 16),
      axis.text      = element_text(size = 16),
      legend.position = "bottom",
      legend.title   = element_text(size = 14),
      legend.text    = element_text(size = 12)
    )
  
  return(PBK_plot)
}

#' plot_PBK_POD 
#'
#' @description Plot a ggplot2  for distribution of BMD, BMDL and BMDU values for the KE from DD genes compared with the PBK prediciton
#'
#' @return Returns ggplot2 object.
#'
#' @import tidyr
#' @import dplyr
#' @import ggplot2
#' @import ggnewscale
#' @import ggrepel
#' @export
#' @noRd
plot_PBK_POD <- function(Enrichment_data_BMD, 
                         output_forward_dosimetry, 
                         exposure_time=40, 
                         exposure_level="X ng/kg/day",
                         y_label = "log10(Concentration (µg/L))",
                         x_label = "Time (days)",
                         free_y_param = F,
                         log_y_scale = TRUE,
                         theme_element_size = 20){
  
  
  Enrichment_data_BMD <- Enrichment_data_BMD |> dplyr::select(Experiment, BMDL, BMD, BMDU, TermID, Ke_description, PFAS, time, Ke_type, level) |> 
    dplyr::distinct() |> 
    dplyr::mutate(BMDL=BMDL |> as.numeric(), BMD=BMD |> as.numeric(), BMDU=BMDU |> as.numeric())

  Enrichment_data_BMD$Experiment = factor(Enrichment_data_BMD$Experiment, levels = Enrichment_data_BMD$Experiment |> unique() |> sort_alphnum())
  Enrichment_data_BMD$Ke_type = factor(Enrichment_data_BMD$Ke_type, levels = c("MolecularInitiatingEvent", "KeyEvent", "AdverseOutcome"))
  Enrichment_data_BMD$level = factor(Enrichment_data_BMD$level, levels = c("Molecular", "Cellular", "Tissue", "Organ", "Individual"))
  Enrichment_data_BMD$time = factor(Enrichment_data_BMD$time, levels = sort(as.numeric(unique(Enrichment_data_BMD$time)), decreasing = F))

  PFOS_mw = 500.13 #g/mol
  PFOA_mw = 414.07 #g/mol
  PFAS_mw = c("PFOS" = PFOS_mw, "PFOA"= PFOA_mw)

  Enrichment_data_BMD <- Enrichment_data_BMD |> dplyr::mutate(BMDL=BMDL*PFAS_mw[PFAS], BMD=BMD*PFAS_mw[PFAS], BMDU=BMDU*PFAS_mw[PFAS])
  KE_long <- Enrichment_data_BMD |> pivot_longer(cols = c(BMDL, BMD, BMDU), names_to = "BMD_type", values_to = "BMD_value")
  KE_long$BMD_type <- factor(KE_long$BMD_type, levels = c("BMDL", "BMD", "BMDU"))

  blue_pal <- c(
    BMDL = "#08306B",  # dark blue
    BMD  = "#4292C6",  # medium blue
    BMDU = "#C6DBEF"   # light blue
  )


  if(log_y_scale){
    KE_long$BMD_value  = log(x = KE_long$BMD_value, base = 10)
  }
  
  p <- ggplot(KE_long, aes(x = time, y = BMD_value, color = BMD_type, fill = BMD_type)) +
    geom_boxplot(
      outlier.shape = NA,
      width = 0.6,
      alpha = 0.7,
      position = position_dodge(width = 0.7)
    ) +
    geom_jitter(
      aes(color = BMD_type),
      alpha = 0.3,
      size = 2,
      position = position_jitterdodge(
        dodge.width = 0.7,
        jitter.width = 0.15
      )
    ) +
    scale_fill_manual(values = blue_pal) +
    scale_color_manual(values = blue_pal) +
    guides(fill = "none") +   # ⬅ remove fill legend only
    facet_wrap(~ PFAS, scales = if (free_y_param) "free_y" else "fixed") +
    theme_minimal(base_size = 20) +
    labs(
      x = x_label,
      y = y_label,
      fill = "KE-level PODs",
      color = "KE-level PODs"
    )

  outFW <- output_forward_dosimetry

  df_PBK <- plot_PBK_simulations(
    outFW,
    scenario_labels = c(exposure_level),
    # chemicals = c("PFOS", "PFOA"),
    chemicals = unique(Enrichment_data_BMD$PFAS),
    return_data = TRUE
  )

  df_PBK_serum <- plot_PBK_simulations_serum(
    outFW,
    scenario_labels = c(exposure_level),
    # chemicals = c("PFOS", "PFOA"),
    chemicals = unique(Enrichment_data_BMD$PFAS),
    return_data = TRUE
  )

  df_PBK_hlines <- df_PBK |> dplyr::filter(Year == exposure_time) |> dplyr::rename(PFAS = Chemical)
  df_PBK_hlines$label <- "Concentration in liver"
  df_PBK_hlines$tissue <- "Liver"

  df_PBK_hlines2 <- df_PBK_serum |> dplyr::filter(Year == 40) |> dplyr::rename(PFAS = Chemical)
  df_PBK_hlines2$label <- "Concentration in Serum"
  df_PBK_hlines2$tissue <- "Serum"

  df_PBK_hlines = rbind(df_PBK_hlines, df_PBK_hlines2)
  df_PBK_hlines$tissue = as.factor(df_PBK_hlines$tissue)

  compartment_pal <- c(
    Liver  = "#1B9E77",  # teal (not green)
    Serum  = "#D8B365"   # sand / tan
  )

  df_PBK_long <- dplyr::bind_rows(
    transform(df_PBK_hlines, y = ConcentrationL, PBK_level = "p5"),
    transform(df_PBK_hlines, y = Concentration,  PBK_level = "p50"),
    transform(df_PBK_hlines, y = ConcentrationU, PBK_level = "p95")
  )
  
  if(log_y_scale){
    df_PBK_long$y  = log(x = df_PBK_long$y, base = 10)
  }
  
  p = p +
    new_scale_color() +
    geom_hline(
      data = df_PBK_long,
      aes(
        yintercept = y,
        color     = tissue,
        linetype  = PBK_level
      ),
      linewidth = 1.1,
      inherit.aes = FALSE
    ) +
    scale_color_manual(name   = "Compartment", values = compartment_pal) +
    scale_linetype_manual(
      name = paste("PBK estimate\n", exposure_level, "\n(40 years)", sep = ""),
      values = c(
        p5  = "dotted",
        p50 = "solid",
        p95 = "dashed"   # .–.–
      ),
      labels = c(
        p5  = "5th quantile",
        p50 = "50th quantile",
        p95 = "95th quantile"
      )
    ) +
    coord_cartesian(clip = "off")

  # if(log_y_scale){
  #   
  #   p <- p + scale_y_log10()
  # }
  
  p  =  p + theme(text = element_text(size = theme_element_size)) 
  
  return(p)
}

