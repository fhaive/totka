
#' dist_plot 
#'
#' @description Generate distribution plot for score
#'
#' @return ggplot object
#'
#' @import dplyr
#' @import ggplot2
#' @export
#' @noRd
dist_plot <- function(KE_annotated, level="all", score="BMD") {
  #KE_annotated <- KE_annotated |> dplyr::select(Experiment, BMD, TermID, Ke_description, PFAS, time, Ke_type, level) |> dplyr::distinct() |> 
  KE_annotated <- KE_annotated |> dplyr::select(Experiment, score={{score}}, TermID, Ke_description, PFAS, time, Ke_type, level) |> dplyr::distinct() |> 
    dplyr::mutate(
      #BMD=BMD |> as.numeric(), 
      #{{score}}:=get(score) |> as.numeric(), 
      score=score |> as.numeric(), 
      Experiment=factor(Experiment, levels=Experiment |> unique() |> sort_alphnum()),
      Ke_type=factor(Ke_type, levels=c("MolecularInitiatingEvent", "KeyEvent", "AdverseOutcome")),
      level=factor(level, levels=c("Molecular", "Cellular", "Tissue", "Organ", "Individual")),
      time=factor(time, levels=time |> unique() |> as.numeric() |> sort() |> as.character())
    )

  PFOS_mw = 500.13 #g/mol 
  PFOA_mw = 414.07 #g/mol 
  PFAS_mw = c("PFOS" = PFOS_mw, "PFOA"= PFOA_mw)

  KE_annotated$score[KE_annotated$PFAS=="PFOS"] = KE_annotated$score[KE_annotated$PFAS=="PFOS"] * PFAS_mw["PFOS"]
  KE_annotated$score[KE_annotated$PFAS=="PFOA"] = KE_annotated$score[KE_annotated$PFAS=="PFOA"] * PFAS_mw["PFOA"]

  if(level=="all"){
    #p.dist = ggplot(KE_annotated, aes(x = time, y = BMD)) +
    p.dist = ggplot(KE_annotated, aes(x = time, y = score)) +
      geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.7, fill = "steelblue") +
      geom_jitter(width = 0.15, alpha = 0.3, size = 2, ) +
      facet_wrap(~PFAS)
  }else{
    #p.dist = ggplot(KE_annotated, aes(x = time, y = BMD, fill = level)) +
    p.dist = ggplot(KE_annotated, aes(x = time, y = score, fill = level)) +
      geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.7) +
      facet_wrap(~PFAS, ncol= 1, scales = "free_y") + 
      scale_fill_brewer(palette = "Blues")
  }
  p.dist <- p.dist +
    theme_minimal(base_size = 20) +
    labs(x = "Time", y = paste0(score, " (µg/L)")) +
    theme(strip.text = element_text(size = 20),)

  return(p.dist)
}

