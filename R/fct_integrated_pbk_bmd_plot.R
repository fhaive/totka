
#' Plot KE activation probabilities over PBK simulation time
#'
#' Creates a tiled heatmap of key-event (KE) activation probabilities across PBK
#' simulation time (x = `Year`, y = `key_event_name`, fill = `Prob`), faceted by
#' `PFAS` and `Day` (columns) and by biological `level` and a PFAS presence label
#' (rows). The PFAS presence label is computed per KE from the occurrence of
#' \code{PFOA} and \code{PFOS} and uses whitespace strings to visually separate
#' facet rows.
#'
#' @details
#' The function merges the input data with the Biological_system_annotations available in the AOPFingerprintR library
#' \code{Biological_system_annotations} using
#' \code{merge(..., by.x = "key_event_name", by.y = "KeyEvent")}.
#' That object contain at least
#' \code{key_event_name} (character) and \code{level} (character) columns.
#'
#' The input \code{df_tissue_specific} must contain (at minimum) the columns:
#' \itemize{
#'   \item \code{KeyEvent} (character) – used for the merge, becomes \code{key_event_name}
#'   \item \code{Year} (numeric/integer) – x-axis
#'   \item \code{PFAS} (character/factor; includes "PFOA" and/or "PFOS") – facet column
#'   \item \code{Prob} (numeric; 0–1) – tile fill
#'   \item \code{Day} (numeric/character) – facet column
#' }
#'
#' The y-axis labels are wrapped using a helper \code{cut_n_words()}. its signature 
#' accept a character vector and an \code{n} argument (number of words).
#'
#' The colour scale is defined via \code{scale_fill_gradientn()} with user-supplied
#' colours, values, breaks, and labels. Faceting uses \code{scales = "free_y"} and
#' \code{space = "free_y"}.
#'
#' @param df_tissue_specific A data frame with the required columns
#'   \code{KeyEvent}, \code{Year}, \code{PFAS}, \code{Prob}, and \code{Day}.
#'   It will be merged with \code{Biological_system_annotations}.
#' @param level_order Character vector giving the desired ordering of the
#'   \code{level} factor (default:
#'   \code{c("Molecular","Cellular","Tissue","Organ","Individual")}).
#' @param pfoa_label Character string used as the PFAS presence label when only
#'   \code{PFOA} is present for a KE (default: \code{""}).
#' @param pfos_label Character string used as the PFAS presence label when only
#'   \code{PFOS} is present for a KE (default: \code{" "} — a single space).
#' @param both_label Character string used as the PFAS presence label when both
#'   \code{PFOA} and \code{PFOS} are present for a KE (default: \code{"   "} — three spaces).
#' @param fill_colours Character vector of colours passed to
#'   \code{ggplot2::scale_fill_gradientn()}.
#' @param prob_values Numeric vector (typically in [0, 1]) passed to
#'   \code{scales::rescale()} to set the colour scale values.
#' @param prob_breaks Numeric vector of breaks for the fill legend.
#' @param prob_labels Character vector of labels corresponding to \code{prob_breaks}.
#' @param x_lab,y_lab Character strings for the x- and y-axis titles.
#' @param labs_fill_text Character string for the fill legend label (via \code{labs(fill = ...)}).
#' @param guide_fill_title Character string for the colourbar title
#'   (via \code{guides(fill = guide_colourbar(title = ...))}).
#' @param y_wrap_n Integer; number of words to keep per line in y-axis labels
#'   (passed to \code{cut_n_words()}).
#' @param axis_text_x_angle,axis_text_x_hjust Numeric; angle and horizontal
#'   justification for x-axis text.
#' @param axis_text_x_size,axis_text_y_size Numeric; font sizes for axis tick labels.
#' @param axis_title_size Numeric; font size for axis titles.
#' @param strip_text_size Numeric; font size for facet strip text.
#' @param legend_title_size,legend_text_size Numeric; font sizes for legend title and text.
#' @param strip_text_y_angle,strip_text_y_hjust Numeric; angle and horizontal
#'   justification for the y-strip text.
#' @param panel_spacing_lines Numeric; vertical panel spacing in "lines"
#'   (passed to \code{grid::unit()}).
#'
#' @return A \code{ggplot} object representing the heatmap.
#'
#' @import dplyr
#' @import tidyr
#' @import AOPfingerprintR
#' @import ggplot2
#' @export
#'
plot_KE_probabilities_over_PBK_time <- function(
    df_tissue_specific,
    # ---- substituted inputs (were hard-coded) ----
    level_order = c("Molecular", "Cellular", "Tissue", "Organ", "Individual"),
    pfoa_label = "",
    pfos_label = " ",
    both_label = "   ",
    fill_colours = c("#BDBDBD","darkblue", "lightblue", "lightgreen", "yellow", "pink", "red"),
    prob_values = c(0,0.1, 0.2, 0.4, 0.6, 0.8, 1),
    prob_breaks = c(0,0.1, 0.2, 0.4, 0.6, 0.8, 1),
    prob_labels = c("0","0.1", "0.2", "0.4", "0.6", "0.8", "1"),
    x_lab = "PBK Years Simulation",
    y_lab = "Key Events",
    labs_fill_text = "BMD",
    guide_fill_title = "Probability of activation",
    y_wrap_n = 10,
    axis_text_x_angle = 45,
    axis_text_x_hjust = 1,
    axis_text_x_size = 12,
    axis_text_y_size = 14,
    axis_title_size  = 14,
    strip_text_size  = 14,
    legend_title_size = 14,
    legend_text_size  = 14,
    strip_text_y_angle = 0,
    strip_text_y_hjust = 0,
    panel_spacing_lines = 1
) {

  Biological_system_annotations = AOPfingerprintR::Biological_system_annotations
  Biological_system_annotations[which(Biological_system_annotations$ke=="Event:164"), "key_event_name"] = "Induction, Liver Dysfunctional Changes by CGA 330050"
  df_tissue_specific$KeyEvent = gsub(pattern = "'",replacement = "",x = df_tissue_specific$KeyEvent)
  
  cut_n_words <- function(x, n = 20) {
    sapply(strsplit(x, "\\s+"), function(words) {
      if (length(words) > n) {
        paste(c(head(words, n), "..."), collapse = " ")
      } else {
        paste(words, collapse = " ")
      }
    }, USE.NAMES = FALSE)
  }
  
  df_tissue_specific = merge(Biological_system_annotations,df_tissue_specific, by.x = "ke", by.y = "Ke", all.y = TRUE)
  df_tissue_specific$key_event_name = gsub(pattern = "'",replacement = "",x = df_tissue_specific$key_event_name)
  
  # Convert 'level' column to a factor with the specified order
  df_tissue_specific$level <- factor(df_tissue_specific$level, levels = level_order, ordered = TRUE)
  # Sort the dataframe by the 'level' column
  df_tissue_specific <- df_tissue_specific %>% arrange(level)
  
  print("str(df_tissue_specific)")
  print(str(df_tissue_specific))
  
  # Step 1: Count occurrences of each KeyEvent per PFAS
  ke_pfas_counts <- df_tissue_specific %>%
    count(key_event_name, PFAS) %>%
    pivot_wider(names_from = PFAS, values_from = n, values_fill = 0)
  
  # Step 2: Create a new column with the label
  # ke_pfas_counts <- ke_pfas_counts %>%
  #   mutate(PFAS_label = case_when(
  #     PFOA > 0 & PFOS == 0 ~ pfoa_label,
  #     PFOS > 0 & PFOA == 0 ~ pfos_label,
  #     PFOA > 0 & PFOS > 0 ~  both_label,
  #     TRUE ~ NA_character_
  #   ))
  
  # Step 3: Join the label back to the original data
  # df_tissue_specific <- df_tissue_specific %>%
  # left_join(ke_pfas_counts %>% dplyr::select(key_event_name, PFAS_label), by = "key_event_name")
  
  gp <- ggplot(df_tissue_specific, aes(x = Year, y = key_event_name, fill = Prob)) +
    geom_tile(color = "white") +
    facet_grid(
      rows = vars(level),#, PFAS_label), # ,level,PFAS_label, organ_tissue.x),
      cols = vars(PFAS, Day),
      scales = "free_y",
      space = "free_y"
    ) +
    # scale_fill_gradient2(low = "darkblue", mid = "orange", high = "red", midpoint = 1) +
    scale_fill_gradientn(
      colours = fill_colours,
      values = prob_values,#scales::rescale(prob_values),
      # breaks = prob_breaks,
      # labels = prob_labels
    ) +
    # wrap by words in the y-axis labels:
    scale_y_discrete(labels = function(labs) cut_n_words(labs, n = y_wrap_n)) +
    theme_minimal() +
    labs(x = x_lab, y = y_lab, fill = labs_fill_text) +
    theme(
      axis.text.x   = element_text(angle = axis_text_x_angle, hjust = axis_text_x_hjust, size = axis_text_x_size),
      axis.text.y   = element_text(size = axis_text_y_size, lineheight = 0.95),
      axis.title.x  = element_text(size = axis_title_size),
      axis.title.y  = element_text(size = axis_title_size),
      strip.text    = element_text(size = strip_text_size),
      legend.title  = element_text(size = legend_title_size),
      legend.text   = element_text(size = legend_text_size),
      strip.text.y  = element_text(angle = strip_text_y_angle, hjust = strip_text_y_hjust, size = strip_text_size),
      panel.spacing.y = unit(panel_spacing_lines, "lines")
    ) +
    guides(fill = guide_colourbar(title = guide_fill_title))
  
  return(gp)
}

