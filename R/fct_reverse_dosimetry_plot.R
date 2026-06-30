
#' plot_reverse_dosimetry_density 
#'
#' @description Plot density for reverse dosimetry
#'
#' @return Returns desnity ggplot.
#'
#' @import dplyr
#' @import ggplot2
#' @import ggridges
#' @import AOPfingerprintR
#' @importFrom reshape2 melt
#' @importFrom stringr str_replace
#' @noRd
plot_reverse_dosimetry_density <- function(reverse_dosimetry_list) {
  RDL <- reverse_dosimetry_list |> reshape2::melt(reverse_dosimetry_list[[1]] |> colnames()) |> dplyr::mutate(PFAS=L1 |> stringr::str_replace(pattern="_.*", replacement="")) |> dplyr::select(-L1)

  Biological_system_annotations <- AOPfingerprintR::Biological_system_annotations |> as.data.frame()
  Biological_system_annotations$organ_tissue[is.na(Biological_system_annotations$organ_tissue)] <- "General"
  rownames(Biological_system_annotations) <- Biological_system_annotations$key_event_name
  
  RDL <- RDL |> mutate(exposure=exposure |> as.character() |> as.numeric() |> suppressWarnings(),
    key_event=key_event |> as.factor() |> suppressWarnings(),
    time=time |> factor(levels=time |> as.numeric() |> unique() |> sort()) |> suppressWarnings()
  )

  plot_df <- RDL
  plot_df <- plot_df |> dplyr::inner_join(Biological_system_annotations, by=c("key_event"="key_event_name")) |> dplyr::mutate(POD=POD |> as.numeric())

  gp <- ggplot(data=plot_df, aes(x = exposure, y = time, fill = key_event)) +
    geom_density_ridges(
      scale = 1.2, rel_min_height = 0.01, alpha = 0.8, color = "grey20"
    ) +
    facet_wrap(~key_event+level, nrow = 4,scales = "free_x") +
    labs(
      x = "Human equivalent exposure (ng/kg/day)", y = "In vitro time (day)"
    ) +
    theme_bw(base_size = 13)

  return(gp)
}

#' plot_reverse_dosimetry_heatmap 
#'
#' @description Plot heatmap for reverse dosimetry
#'
#' @return Returns heatmap ggplot.
#'
#' @import dplyr
#' @import tidyr
#' @import reshape2
#' @importFrom stringr str_replace
#' @export
#' @noRd
plot_reverse_dosimetry_heatmap <- function(reverse_dosimetry_list_summarized, score_col="exposure", name="Exposure ng/kg/day", row_title="Key Events", column_title="Experiment", col_fun=NULL, row_ha=NULL, col_ha=NULL, row_names_gp=NULL, column_names_gp=NULL, legend.padding=unit(c(2, 2, 2, 100),"mm"), toDraw=TRUE) {
  RD <- reverse_dosimetry_list_summarized |> reshape2::melt(id.vars=reverse_dosimetry_list_summarized[[1]] |> colnames()) |> dplyr::mutate(Experiment=L1 |> stringr::str_replace(pattern="_[:alpha:]*$", replacement="")) |> dplyr::mutate(Experiment=factor(Experiment, levels=Experiment |> sort_alphnum()))

  Biological_system_annotations <- AOPfingerprintR::Biological_system_annotations |> as.data.frame()
  Biological_system_annotations$organ_tissue[is.na(Biological_system_annotations$organ_tissue)] <- "General"
  rownames(Biological_system_annotations) <- Biological_system_annotations$key_event_name

  RD <- RD |> dplyr::inner_join(Biological_system_annotations, by=c("key_event"="key_event_name"))
  RD <- RD |> dplyr::mutate(level=level |> factor(levels=c("Molecular","Cellular","Tissue","Organ","Individual")), pfas=Experiment |> as.vector() |> stringr::str_replace(pattern="_.*", replacement=""), time=time |> as.numeric()) 

  mat <- reshape2::acast(RD, key_event~Experiment, value.var = {{score_col}})
  mat[is.na(mat)] <- 0
  colnames_mat <- colnames(mat)
  pfas_type <- sub("_.*", "", colnames_mat)              # "PFOS" / "PFOA"
  time_point <- sub(".*_", "", colnames_mat)             # keep as character for discrete

  if(ncol(mat)>1){
    mat <- mat[,colnames_mat |> sort_alphnum()]
  }

  ## Define discrete colors for time points
  #time_levels <- unique(time_point)
  #time_colors <- setNames(
  #  RColorBrewer::brewer.pal(n = length(time_levels), name = "Set2"), 
  #  time_levels
  #)

  ## Define column annotation
  #col_ha <- HeatmapAnnotation(
  #  PFAS=pfas_type,
  #  Time=factor(time_point, levels = time_levels |> as.numeric() |> sort()),
  #  col=list(
  #    PFAS=c("PFOS"="#1f78b4", "PFOA"="#33a02c"),  # categorical PFAS colors
  #    Time=time_colors                                 # categorical time colors
  #  )
  #)

  # pick distinct colors for each organ/tissue
  level_cols <- setNames(
    c("#1b9e77","#d95f02","#7570b3","#e7298a","#08306b"),
    c("Organ","Cellular","Molecular","Individual","Tissue")
  )

  row_ha <- rowAnnotation(
    Level=Biological_system_annotations[rownames(mat),"level"],
    col=list(Level=level_cols)
    # Organ = Biological_system_annotations[rownames(mat),"organ_tissue"]
  )

  #rng <- range(mat, na.rm = TRUE)

  #if (rng[1] < 0 && rng[2] > 0) {
  #  # mixed negatives & positives: diverging with 0 = white
  #  col_fun <- circlize::colorRamp2(
  #    c(rng[1], 0, rng[2]),
  #    c("#2c7bb6", "white", "#d7191c")
  #  )
  #} else if (rng[1] >= 0) {
  #  # all non-negative
  #  rng <- range(mat[mat>0], na.rm = TRUE)
  #  col_fun <- circlize::colorRamp2(
  #    c(0, rng[1], rng[2]),  # breakpoints: 0, low, max
  #    c("white","#FF7C7C", "#2c7bb6") # white, light blue, dark blue
  #  )
  #} else {
  #  # all non-positive: blue at min up to white at 0
  #  col_fun <- circlize::colorRamp2(
  #    c(rng[1], 0),
  #    c("#2c7bb6", "white")
  #  )
  #}

  #ht <- Heatmap(
  #  mat,
  #  name = "Exposure ng/kg/day",
  #  col = col_fun,
  #  na_col = "grey90",
  #  # row_split = organ,
  #  left_annotation = row_ha,
  #  top_annotation = col_ha,
  #  cluster_rows = TRUE,
  #  cluster_columns = F,
  #  show_row_names = TRUE,
  #  show_column_names = TRUE,
  #  show_row_dend = FALSE,
  #  show_column_dend = FALSE,
  #  row_title = "Key Events",
  #  column_title = "Experiment"
  #)

  #ht <- draw(
  #  ht,
  #  heatmap_legend_side = "left",
  #  annotation_legend_side = "top",
  #  #padding = unit(c(2, 2, 2, 290), "mm")
  #  padding = unit(c(2, 2, 2, 100), "mm")
  #)

  ht <- plot_heatmap(mat=mat, name=name, row_title=row_title, column_title=column_title, col_fun=col_fun, row_ha=row_ha, col_ha=col_ha, row_names_gp=row_names_gp, column_names_gp=column_names_gp, legend.padding=legend.padding, toDraw=toDraw)
  return(ht)
}

