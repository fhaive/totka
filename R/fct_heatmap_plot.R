
#' plot_heatmap 
#'
#' @description A function to plot ComplexHeatmap heatmap
#'
#' @return ComplexHeatmap object
#'
#' @import tidyr
#' @import ComplexHeatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar
#'
#' @export
#' @noRd
plot_heatmap <- function(mat, pfas_type=mat|>colnames()|>stringr::str_replace(pattern="_.*", replacement=""), time_point=mat|>colnames()|>stringr::str_replace(pattern=".*_", replacement=""), name="BMD", row_title="Genes", column_title="Experiment", col_fun=NULL, row_ha=NULL, col_ha=NULL, row_names_gp=NULL, column_names_gp=NULL, legend.padding=unit(c(2, 2, 2, 20),"mm"), toDraw=TRUE) {
  mat[is.na(mat)] = 0
  colnames_mat <- colnames(mat)
  if(ncol(mat)>1){
    mat <- mat[,colnames_mat |> sort_alphnum()]
  }

  # Define discrete colors for time points
  time_levels <- unique(time_point) |> as.numeric() |> sort() |> as.character()

  if(length(time_levels)<3){
    col.vec <- RColorBrewer::brewer.pal(n = 3, name = "Set2")
    col.vec <-  col.vec[1:length(time_levels)]
  }else{
    col.vec <- RColorBrewer::brewer.pal(n = length(time_levels), name = "Set2")
  }

  time_colors <- setNames(
    #RColorBrewer::brewer.pal(n = length(time_levels), name = "Set2"), 
    col.vec, 
    time_levels
  )

  # Define column annotation
  if(is.null(col_ha)){
    col_ha <- ComplexHeatmap::HeatmapAnnotation(
      PFAS = pfas_type,
      Time=factor(time_point, levels = time_levels),
      col = list(
        PFAS = c("PFOS" = "#1f78b4", "PFOA" = "#33a02c"),  # categorical PFAS colors
        Time = time_colors                                 # categorical time colors
      ),
      annotation_name_gp = grid::gpar(fontsize = 18)
    )
  }
  
  if(is.null(col_fun)){
    rng <- range(mat, na.rm = TRUE)

    if (rng[1] < 0 && rng[2] > 0) {
      # mixed negatives & positives: diverging with 0 = white
      col_fun <- circlize::colorRamp2(
        c(rng[1], 0, rng[2]),
        c("#2c7bb6", "white", "#d7191c")
      )
    } else if (rng[1] >= 0) {
      # all non-negative
      rng <- range(mat[mat>0], na.rm = TRUE)
      col_fun <- circlize::colorRamp2(
        c(0, rng[1], rng[2]),  # breakpoints: 0, low, max
        c("white","#FF7C7C", "#2c7bb6") # white, light blue, dark blue
      )
    } else {
      # all non-positive: blue at min up to white at 0
      col_fun <- circlize::colorRamp2(
        c(rng[1], 0),
        c("#2c7bb6", "white")
      )
    }
  }

  if(is.null(row_names_gp)){
    row_names_gp = grid::gpar(fontsize = 18)
  }

  if(is.null(column_names_gp)){
    column_names_gp = grid::gpar(fontsize = 18)
  }
  
  ht <- ComplexHeatmap::Heatmap(
    mat,
    name = name,
    col = col_fun,
    na_col = "grey90",
    left_annotation = row_ha,
    top_annotation = col_ha,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = row_names_gp,
    column_names_gp = column_names_gp,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    row_title = row_title,
    column_title = column_title
  )

  if(toDraw){
    ht <- ComplexHeatmap::draw(
      ht,
      heatmap_legend_side = "left",
      annotation_legend_side = "top",
      padding = legend.padding
    )
  }

  return(ht)
}

#' plot_heatmap_preprocess 
#'
#' @description A function to plot ComplexHeatmap heatmap for preprocessed input data Differential Expression analysis OR BMDx results table
#'
#' @return ComplexHeatmap object
#'
#' @import dplyr
#' @importFrom reshape2 acast
#'
#' @export
#' @noRd
plot_heatmap_preprocess <- function(DF, score_col="log2FoldChange", name="log2FoldChange", row_title="Genes", column_title="Experiment", col_fun=NULL, row_ha=NULL, col_ha=NULL, row_names_gp=NULL, column_names_gp=NULL, legend.padding=unit(c(2, 2, 2, 20),"mm"), toDraw=TRUE) {
  DF <- DF |> dplyr::mutate(Experiment=paste0(Experiment, "_", time)) |> dplyr::select(Experiment, Feature, {{score_col}}) |> dplyr::distinct()
  mat <- DF |> reshape2::acast(formula = Feature~Experiment, value.var={{score_col}})
  mat[is.na(mat)] <- 0
  ht <- plot_heatmap(mat=mat, name=name, row_title=row_title, column_title=column_title, col_fun=col_fun, row_ha=row_ha, col_ha=col_ha, row_names_gp=row_names_gp, column_names_gp=column_names_gp, legend.padding=legend.padding, toDraw=toDraw)
  return(ht)
}

#' plot_heatmap_ke 
#'
#' @description A function to plot ComplexHeatmap heatmap for key events
#'
#' @return ComplexHeatmap object
#'
#' @import dplyr
#' @import tidyr
#' @importFrom randomcoloR distinctColorPalette
#' @importFrom colorspace darken
#'
#' @export
#' @noRd
plot_heatmap_ke <- function(DF, score_col="log2FoldChange", name="log2FoldChange", row_title="KE (split & annotated by organ/tissue)", column_title="Experiment", legend.padding=unit(c(2, 2, 2, 290),"mm"), toDraw=TRUE) {
  DF <- DF |> 
    dplyr::select(TermID, Experiment, organ_tissue, {{score_col}}, Ke_description) |> 
    dplyr::mutate(organ_tissue=ifelse(is.na(organ_tissue),"uncategorized",organ_tissue), {{score_col}}:=get(score_col) |> as.numeric()) |> 
    dplyr::distinct() |> 
    dplyr::arrange(get(score_col)) |> 
    dplyr::mutate(Ke_description=factor(Ke_description, levels=Ke_description |> unique()))

  df_wide <- DF %>%
    dplyr::select(organ_tissue, Ke_description, Experiment, {{score_col}}) %>%
    group_by(organ_tissue, Ke_description, Experiment) %>%
    summarise({{score_col}}:=mean(get(score_col), na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Experiment, values_from = {{score_col}})

  # Extract the matrix and the row split vector
  mat <- df_wide %>%
    dplyr::select(-organ_tissue, -Ke_description) %>%
    as.data.frame() %>%
    as.matrix()

  # Row names = KE descriptions (ensure uniqueness if same KE appears in >1 organ)
  rownames(mat) <- make.unique(as.character(df_wide$Ke_description))

  # vector aligned to rows of `mat`
  organ <- factor(df_wide$organ_tissue)
  organ_lvls <- levels(organ)

  # pick distinct colors for each organ/tissue
  organ_cols <- setNames(
    randomcoloR::distinctColorPalette(length(organ_lvls)) |> colorspace::darken(amount=0.3),
    organ_lvls
  )

  # row annotation (a colored strip + legend)
  row_ha <- rowAnnotation(
    organ_tissue = organ,
    col = list(organ_tissue = organ_cols),
    annotation_legend_param = list(title = "Organ/Tissue"),
    annotation_name_side = NULL
  )

  row_names_gp = gpar(col = organ_cols[organ])   # color KE labels by organ

  mat[is.na(mat)] <- 0
  ht <- plot_heatmap(mat=mat, name=name, row_title=row_title, column_title=column_title, row_ha=row_ha, row_names_gp=row_names_gp, legend.padding=legend.padding, toDraw=toDraw)
  return(ht)
}

