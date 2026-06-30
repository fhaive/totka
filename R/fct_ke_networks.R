
### ----------------------------------------------------------
### SIMPLE PATH SUPPORT SCORING FUNCTION
### ----------------------------------------------------------
#' score_path
#'
#' @description Simple path support scoring function
#'
#' @return Returns list of normalized scores.
#'
score_path <- function(path_nodes,
                       DD, 
                       DEG,
                       BMD,
                       PBK, 
                       DEG_score, 
                       Binding,
                       w_BMD = 2, 
                       w_PBK = 1.5, 
                       w_DEG = 1, 
                       w_Binding = 0.5,
                       alpha = 0.2,
                       normalize = NULL) {

  w_PBK = 0.8 #0.6
  w_DEG = 0.3
  w_Binding = 0.2
  alpha = 0.6 #0.8

  # Subset evidence to only the nodes in the path
  dd  <- DD[path_nodes]
  deg <- DEG[path_nodes]
  pbk <- PBK[path_nodes]
  dar <- DEG_score[path_nodes]
  bind <- Binding[path_nodes]
  
  # Replace missing values with 0
  if(!is.null(dd)){
    dd[is.na(dd)] <- 0
  }else{
    dd = 0
  }
  
  if(!is.null(deg)){
    deg[is.na(deg)] <- 0
  }else{
    deg = 0
  }
  
  if(!is.null(pbk)){
    pbk[is.na(pbk)] <- 0
  }else{
    pbk = 0
  }
  
  if(!is.null(dar)){
    dar[is.na(dar)] <- 0
  }else{
    dar = 0
  }
  
  if(!is.null(bind)){
    bind[is.na(bind)] <- 0
  }else{
    bind = 0
  }
  
  # Normalize each evidence type within the path
  dd   <- if (max(dd)==0)   dd   else dd   /  max(dd)
  deg  <- if (max(deg)==0)  deg  else deg  /  max(deg)
  pbk  <- if (max(pbk)==0)  pbk  else pbk  /  max(pbk)
  dar  <- if (max(dar)==0)  dar  else dar  /  max(dar)
  bind <- if (max(bind)==0) bind else bind /  max(bind)
  
  # Node-level evidence contributions
  S_pbek   <- w_PBK    * pbk  * dd
  S_deg  <- w_DEG  * dar  * deg
  S_bind   <- w_Binding * bind
  
  # Total node score
  S_nodes <-  S_pbek + S_deg + S_bind
  names(S_nodes) = path_nodes
  
  # Raw path score
  S_raw <- sum(S_nodes)
  
  # Coverage: fraction of nodes with at least one evidence
  Ecount <- (pbk > 0) + (dar > 0) + (bind > 0)
  names(Ecount) = path_nodes
  coverage <- sum(Ecount > 0) / length(Ecount)
  
  # Multi-evidence bonus
  multi_bonus <- mean(pmax(Ecount - 1, 0))
  
  if(is.null(normalize)){
    denom = 1
  }else{
    if(normalize == "length"){
      denom <- length(path_nodes) # normalize by path length
    }
    if(normalize == "support"){
      supported <- sum(Ecount > 0) # normalize by nr of nodes with evidences
      denom <- if (supported > 0) supported else 1
    }
  }
 
  # Final score
  S_final <- (S_raw/denom) * coverage * (1 + alpha * multi_bonus)
  
  # if(!is.null(BMD)){
  #   S_final = S_final * (1-BMD)
  # }
  
  list(
    dd = dd,
    deg = deg,
    pbk=pbk,
    dar= dar,
    bind = bind,
    path_score = S_final,
    raw = S_raw,
    coverage = coverage,
    Ecount = Ecount,
    multi_bonus = multi_bonus,
    node_scores = setNames(S_nodes, path_nodes)
  )
}


# -----------------------------
# Path computation (unchanged logic, including the original pp/mx behavior)
# -----------------------------
#' compute_paths_for_experiment 
#'
#' @description Path computation (unchanged logic, including the original pp/mx behavior) 
#'
#' @return Returns list of paths to target node in network.
#'
#' @import dplyr
#' @import igraph
compute_paths_for_experiment <- function(
    expi, g_main, ke_list_exp, docking_list, Enrichment_data_BMD, ke_all_deg_list,
    target_node = "Increased, Liver Steatosis",
    normalize = NULL,
    PBK_time = 40
) {
  in_degree <- igraph::degree(graph = g_main, mode = "in")
  zero_in_degree <- in_degree[in_degree == 0]
  from <- zero_in_degree
  
  is_max <- 0
  opt_path <- NULL
  path_list <- list()
  
  for (i in 1:length(from)) {
    paths <- igraph::all_simple_paths(g_main, from = names(from)[i], to = target_node, mode = "out")
    
    if (length(paths) > 0) {
      # path_sum <- c()
      # path_sum_perc <- c()
      # path_l <- c()
      
      path_scores = c()
      path_statistic_list = c()
      
      for (pp in 1:length(paths)) {
        path_nodes = names(paths[[pp]])
        
        # KE with dose dependent gene mapping in the path
        Enrichment_data_BMD_df <- unique(Enrichment_data_BMD[, c("Ke_description", "Experiment", "relevantGenesInGeneSet","BMD")])
        idx <- Enrichment_data_BMD_df$Ke_description %in% path_nodes &Enrichment_data_BMD_df$Experiment %in% expi
        n_DD_genes_in_enriched_KD <- Enrichment_data_BMD_df$relevantGenesInGeneSet[idx]
        names(n_DD_genes_in_enriched_KD) <- Enrichment_data_BMD_df$Ke_description[idx]
        
        BMD_genes_in_enriched_KD <- as.numeric(Enrichment_data_BMD_df$BMD[idx])
        names(BMD_genes_in_enriched_KD) <- Enrichment_data_BMD_df$Ke_description[idx]
        
        # if(length(BMD_genes_in_enriched_KD)>0){
        #   max_BMD_across_KE = max(BMD_genes_in_enriched_KD)
        #   F <- ecdf(as.numeric(Enrichment_data_BMD_df$BMD))
        #   quantile_BMD <- F(max_BMD_across_KE) # the strongest effect within the group
        # }else{
        #   quantile_BMD = NULL
        # }
        
        # KE with average probability of activation across the PBK years of simulatom
        # max_pbk_omics_prob_in_path <- unlist(lapply(ke_list_exp[[expi]][names(ke_list_exp[[expi]]) %in% path_nodes], median))
        
        #PBK probability at a specific PBKsimulation year for which the network is built
        max_pbk_omics_prob_in_path <- unlist(lapply(ke_list_exp[[expi]][names(ke_list_exp[[expi]]) %in% path_nodes], function(probs_over_years) as.numeric(probs_over_years[PBK_time])))
        
        # KE with score coming from differential analysis -log(pvalue) x logFC
        idx <- ke_all_deg_list[[expi]]$Ke_description %in%path_nodes
        deg_score_in_path <- ke_all_deg_list[[expi]]$prob[idx]
        names(deg_score_in_path) <- ke_all_deg_list[[expi]]$Ke_description[idx]
        
        # Nr of differentially expressed genes in KE
        idx <- ke_all_deg_list[[expi]]$Ke_description %in%path_nodes
        n_DEG_genes_in_enriched_KD <- ke_all_deg_list[[expi]]$relevantGenesInGeneSet[idx]
        names(n_DEG_genes_in_enriched_KD) <- ke_all_deg_list[[expi]]$Ke_description[idx]
        
        # idx <- docking_list[[expi]]$Ke_description %in% path_nodes
        # max_binding_score_in_path <- abs(docking_list[[expi]]$Binding[idx])
        # names(max_binding_score_in_path) <- docking_list[[expi]]$Ke_description[idx]
        
        idx <- docking_list[[expi]]$ke %in% path_nodes
        df = as.data.frame(docking_list[[expi]][idx,])
        
        # columns in dataframe df PFAS, KE, KE_description, GENE, Binding
        if(nrow(df)>0){
          res_by_ke <- df %>%
            group_by(ke) %>%
            # count how many genes (rows) map to this KE
            mutate(gene_count = n_distinct(GENE)) %>%  # or use n() if duplicates are meaningful
            # keep the row(s) with the minimum Binding per KE_description
            filter(Binding == min(Binding, na.rm = TRUE)) %>%
            # if there are ties (same min Binding), keep just one deterministically
            slice(1L) %>%
            ungroup()
          
          max_binding_score_in_path <- abs(res_by_ke$Binding)
          names(max_binding_score_in_path) <-res_by_ke$ke
        }else{
          max_binding_score_in_path <- NULL
        }
        
        # (use your actual vectors)
        res <- score_path(path_nodes,
                          DD = n_DD_genes_in_enriched_KD,
                          DEG = n_DEG_genes_in_enriched_KD,
                          BMD = NULL,
                          PBK = max_pbk_omics_prob_in_path,
                          DEG_score = deg_score_in_path,
                          Binding = max_binding_score_in_path,
                          normalize = normalize)
        
        path_scores = c(path_scores, res$path_score)
        path_statistic_list[[pp]] = res
     
      }
      
      current <- path_scores
      mx <- which.max(current)
      
      # Preserve original behavior (uses `pp` after loop — as in your code)
      if (current[mx] > is_max) {
        is_max <- current[mx]
        opt_path <- names(paths[[pp]])
      }
      
      path_list[[names(from)[i]]] <- list(
        "paths" = paths,
        path_scores = path_scores,
        path_statistic_list = path_statistic_list,
        mx = mx
      )
    }
  }
  
  list(path_list = path_list, opt_path = opt_path, is_max = is_max)
}


#' enlarge_ke_set 
#'
#' @description Enlarge an input set of Key Events by finding connected AO and MIE Events 
#'
#' @return Returns enlarged Key Events set.
#'
#' @import igraph
#' @import AOPfingerprintR
#' @export
#' @noRd
enlarge_ke_set <- function(ke_set, max_path_length=20, n_AOs=10, n_MIEs=10, mode="out"){
  KEKE_net <- igraph::upgrade_graph(AOPfingerprintR::KEKE_net)
  mask <- !is.na(E(KEKE_net)$type.r) & trimws(E(KEKE_net)$type.r) == "DIRECTLY_LEADS_TO"
  KEKE_net <- subgraph_from_edges(KEKE_net, eids = E(KEKE_net)[mask], delete.vertices = FALSE)
  MIE <- AOPfingerprintR::aop_ke_table_hure |> dplyr::filter(Ke_type=="MolecularInitiatingEvent") |> dplyr::pull(Ke) |> unique()
  AO <- AOPfingerprintR::aop_ke_table_hure |> dplyr::filter(Ke_type=="AdverseOutcome") |> dplyr::pull(Ke) |> unique()
  closest_AOs <- AOPfingerprintR::compute_the_closest_AOs(ke_set,
    KEKE_net, AO = AO, max_path_length = max_path_length,
    n_AOs = n_AOs, mode = mode)
  closest_MIEs <- AOPfingerprintR::compute_the_closest_MIEs(ke_set,
    KEKE_net, MIE = MIE, max_path_length = max_path_length,
    n_MIEs = n_MIEs, mode = mode)
  events <- c(ke_set, closest_AOs, closest_MIEs) |> unique()
  return(events)
}

#' get_igraph 
#'
#' @description Given a set of Key Events infer a sub graph from AOPfingerprintR::KEKE_net
#'
#' @return Returns inferred igraph for Key Events set.
#'
#' @import igraph
#' @import AOPfingerprintR
#' @export
#' @noRd
get_igraph <- function(ke_set, enlarge_ke_selection=FALSE, max_path_length=20, n_AOs=10, n_MIEs=10, mode="out"){
    KEKE_net <- igraph::upgrade_graph(AOPfingerprintR::KEKE_net)
    mask <- !is.na(E(KEKE_net)$type.r) & trimws(E(KEKE_net)$type.r) == "DIRECTLY_LEADS_TO"
    KEKE_net <- subgraph_from_edges(KEKE_net, eids = E(KEKE_net)[mask], delete.vertices = FALSE)
    
    interesting_vertex <- ke_set |> unique()
    if (enlarge_ke_selection) {
        events <- enlarge_ke_set(ke_set=interesting_vertex, max_path_length=max_path_length, n_AOs=n_AOs, n_MIEs=n_MIEs, mode=mode)
    }else {
        events <- interesting_vertex
    }
    events <- events[!is.na(events)]
    events <- events[events %in% AOPfingerprintR::aop_ke_table_hure$Ke]
    selegoG <- igraph::induced_subgraph(KEKE_net, igraph::V(KEKE_net)[which(igraph::V(KEKE_net)$name %in% events)])
    return(selegoG)
}

#' get_visnet 
#'
#' @description Given a list of node and edge data frames, generate a vizNetwork
#'
#' @return Returns visNetwork plot.
#'
#' @import dplyr
#' @import igraph
#' @import visNetwork
#' @export
#' @noRd
get_visnet <- function(nodes, edges, group_by="aop", col_pal=NULL){
  # --- Legend: border color (evidence) ---
  if("evidence" %in% colnames(nodes)) {
    lnodes_border <- nodes |>
      dplyr::select(evidence, color.border) |>
      dplyr::distinct() |>
      dplyr::arrange(evidence) |>
      dplyr::mutate(
        label = paste0("Evidences: ", evidence),
        shape = "dot",
        title = "Border color",
        color.background = "white",   # essential so legend item is visible
        borderWidth = ifelse(evidence==0, 1, 6)   # thicker legend borders
      )
    
    no_prob_color_border = lnodes_border$color.border[lnodes_border$evidence==0]
  }

  # ---- Build background-color legend: LOW / MID / HIGH PBK_prob ----
  
  # extract valid PBK values
  pbk_vals <- nodes$PBK_prob[!is.na(nodes$PBK_prob)]
  
  
  # compute representative values
  
  zero_prob_color <- "#BFBFBF"
  
  pbk_min <- 1
  pbk_med <- 50
  pbk_max <- 100
  
  no_prob_color_bg = nodes |> dplyr::filter(is.na(PBK_prob)) |> dplyr::pull(color.background) |> unique()

  # reference table with labels
  lnodes_bg <- data.frame(
    label = c("No Probability","Zero Probability","Low Probability", "Mid Probability", "High Probability"),
    color.background = c(no_prob_color_bg,
                         
                         zero_prob_color,
                         
                         col_pal(pbk_min),
                         col_pal(pbk_med),
                         col_pal(pbk_max)),
    color.border = c("#f2f2f2",
                     
                     zero_prob_color,
                     
                         col_pal(pbk_min),
                         col_pal(pbk_med),
                         col_pal(pbk_max)),
    borderWidth = c(0.2,0,0,0,0),
    shadow = TRUE,
    shape = "dot"
  )
  
  # --- Legend: node-type (shapes) ---
  lnodes_type <- data.frame(
    label = c("MolecularInitiatingEvent", "KeyEvent", "AdverseOutcome"),
    shape = c("triangle", "square", "star"),
    color.background = "white",
    color.border = "black",
    title = "Node types"
  )
  
  # --- Combine all legend entries ---
  if("evidence" %in% colnames(nodes)) {
    legend_nodes <- dplyr::bind_rows(lnodes_type, lnodes_border, lnodes_bg)
  } else{
    legend_nodes <- dplyr::bind_rows(lnodes_type, lnodes_bg)
  }

  # --- Build network ---
  vn <- visNetwork::visNetwork(nodes, edges) %>%
    visNetwork::visLegend(
      position = "right",
      addNodes = legend_nodes,
      useGroups = FALSE
    ) %>%
    visNetwork::visOptions(
      selectedBy = list(variable = group_by, multiple = TRUE),
      highlightNearest = FALSE
    ) %>%
    visNetwork::visNodes(shadow = FALSE) %>%
    visNetwork::visEdges(arrows = 'to', color = "black")

  return(vn)
}

#' get_visnet_reverse_dosimetry
#'
#' @description Given a list of node and edge data frames, generate a vizNetwork
#'
#' @return Returns visNetwork plot.
#'
#' @import grid
#' @import visNetwork
#' @import ggplot2
#' @importFrom cowplot get_legend
#' @export
#' @noRd
get_visnet_reverse_dosimetry <- function(nodes, edges, group_by="aop", col_pal=NULL){
  # ---- Build background-color legend: LOW / MID / HIGH exposure ----
  
  # extract valid exposure values
  score_vals <- nodes$exposure[!is.na(nodes$exposure)]
  
  # compute representative values
  score_min <- min(score_vals)
  medIdx <- (score_vals - median(score_vals)) |> abs() |> which.min()
  score_med <- score_vals[medIdx]
  score_max <- max(score_vals)
  
  # gradient legend
  x <- LETTERS[1:20]
  y <- paste0("var", seq(1,20))
  data <- expand.grid(X=x, Y=y)
  data$Exposure <- runif(400, score_min, score_max)
  breaks.vec <- quantile(data$Exposure, probs=c(0,0.25,0.5,0.75,1)) |> round() |> as.vector()

  # Heatmap
  gplot <- ggplot(data, aes(X, Y, fill= Exposure)) +
    geom_tile() +
    scale_fill_gradient(low=col_pal(score_min), high=col_pal(score_max), breaks=breaks.vec)

  gplot.legend <- cowplot::get_legend(gplot)

  # temp gradient image file
  #tempFile <- tempfile(pattern="gradient_legend_", tmpdir="inst/app/www", fileext=".tiff")
  tempFile <- tempfile(pattern="gradient_legend_", tmpdir=tempdir(), fileext=".tiff")
  #tempFile <- tempfile(pattern="gradient_legend_", fileext=".tiff")

  # Create new plot window
  # Draw Only legend
  tiff(filename=tempFile, res=600, height=1000)
  grid::grid.newpage()
  grid::grid.draw(gplot.legend)
  dev.off()

  if(exists("r") && is.list(r)){
    if(is.null(r$tempFiles)){
      r$tempFiles <- c(tempFile) 
    }else{
      r$tempFiles <- r$tempFiles |> append(tempFile) 
    }
  }

  ##gradientFile <- file.path("www", tempFile |> basename())
  gradientFile <- file.path("tmpDIR", tempFile |> basename())
  #gradientFile <- tempFile
  lnodes_gradient <- data.frame(
    label = "",
    shape = "image",
    color.background = "",
    color.border = "",
    title = "",
    shadow = FALSE,
    image = gradientFile,
    size = 40
  )

  # --- Legend: node-type (shapes) ---
  lnodes_type <- data.frame(
    label = c("MolecularInitiatingEvent", "KeyEvent", "AdverseOutcome"),
    shape = c("triangle", "square", "star"),
    color.background = "white",
    color.border = "black",
    title = "Node types",
    shadow = FALSE,
    image = "",
    size = 15
  )

  # --- Legend: node-type (shapes) ---
  lnodes_gap <- data.frame(
    label = "",
    shape = "dot",
    color.background = "white",
    color.border = "white",
    title = "",
    shadow = FALSE,
    image = "",
    size = 15
  )
  
  # --- Combine all legend entries ---
  legend_nodes <- dplyr::bind_rows(
    lnodes_type,
    lnodes_gap,
    lnodes_gradient)
    #lnodes_bg)
  #print("legend_nodes")
  #print(legend_nodes)

  # --- Build network ---
  vn <- visNetwork::visNetwork(nodes, edges) %>%
    visNetwork::visLegend(
      position = "right",
      addNodes = legend_nodes,
      useGroups = FALSE
    ) %>%
    visNetwork::visOptions(
      selectedBy = list(variable = group_by, multiple = TRUE),
      highlightNearest = FALSE
    ) %>%
    visNetwork::visNodes(shadow = FALSE) %>%
    visNetwork::visEdges(arrows = 'to', color = "black")

  return(vn)
}


#' get_mean_pbk 
#'
#' @description Get KE time series, Time Point(i) table of PBK probability (mean over samples)
#'
#' @return Returns a list of KE time series and mean probaility table for each time point.
#'
#' @import dplyr
#' @noRd
get_mean_pbk <- function(mi3d){
  ## Extract 3D matrix: KEs x samples x timepoints
  #if (is.null(mi3d)) {
  #  warning("No matrix found for experiment '", expi, "'. Skipping.")
  #  next
  #}
  if (length(dim(mi3d)) != 3) {
    stop("Matrix for '", expi, "' must be 3-dimensional (KEs x samples x timepoints).")
  }

  n_time <- dim(mi3d)[3]
  # Build mean-over-samples for each timepoint (KEs x timepoints)
  mean_over_samples <- NULL
  for (tp in seq_len(n_time)) {
    # Row means across samples for timepoint tp
    rms <- rowMeans(mi3d[, , tp], na.rm = TRUE)
    mean_over_samples <- cbind(mean_over_samples, rms)
  }
  mi <- mean_over_samples
  rownames(mi) <- rownames(mi3d) # ensure KE names present
  # set colnames
  colnames(mi) <- seq_len(n_time)

  # Per-KE time series list
  ke_list <- vector(mode = "list", length = nrow(mi))
  names(ke_list) <- rownames(mi)
  for (i in seq_len(nrow(mi))) {
    ke_list[[rownames(mi)[i]]] <- mi[i, ]
  }
  
  return(list(mi=mi, ke_list=ke_list))
}


#' get_pbk_df 
#'
#' @description Get Time Point(i) specific table of PBK probability from 2D probability matrix
#'
#' @return Returns probaility table for specific timepoint.
#'
#' @import dplyr
#' @noRd
get_pbk_df <- function(mi, PBK_time){
  # Base DF with prob at PBK_time
  if (PBK_time < 1 || PBK_time > ncol(mi)) {
    stop("PBK_time=", PBK_time, " is out of range. Available timepoints: 1..", ncol(mi))
  }
  df <- data.frame(
    ke   = rownames(mi),
    PBK_prob = mi[, PBK_time],
    stringsAsFactors = FALSE
  )
  return(df)
}

#' get_deg_list 
#'
#' @description Get formatted list of data.frames for each experiment from Expression KE Enrichment Data
#'
#' @return Returns a list formatted data.frames.
#'
#' @import dplyr
#' @import purrr
#' @noRd
get_deg_list <- function(Enrichment_data_DEG){
  expi_vec <- Enrichment_data_DEG |> dplyr::pull(expi) |> unique()

  # DEG enrichment: aggregate by KE (min padj)
  df_deg_list <- expi_vec |> purrr::map(\(x){
    Enrichment_data_DEG |> dplyr::filter(expi=={{x}}) |> dplyr::select(Ke_description, relevantGenesInGeneSet, padj, log2FoldChange) |> 
      dplyr::group_by(Ke_description) |>
      dplyr::slice_min(order_by=padj, n=1, with_ties=FALSE) |> 
      dplyr::ungroup() |>
      as.data.frame()
  }) |> setNames(expi_vec)

 df_deg_list <- df_deg_list |> purrr::map(\(x){
    x |> dplyr::mutate(pi_value=-log10(as.numeric(padj)) * abs(as.numeric(log2FoldChange)), prob=pi_value)
  })

  ke_all_deg_list <- df_deg_list

  # NOTE: reuse 'prob' column for padj to align with downstream merge
  df_deg_list <- df_deg_list |> purrr::map(\(x){
    x |> dplyr::select(ke=Ke_description, pi_value, DE_padj=padj, log2FoldChange) |> 
      dplyr::mutate(log2FoldChange=log2FoldChange |> as.numeric()) |> 
      dplyr::distinct()
  })
  
  return(list(df_deg_list=df_deg_list, ke_all_deg_list=ke_all_deg_list))
}

#' get_bmd_list 
#'
#' @description Get formatted list of data.frames for each experiment from BMD KE Enrichment Data
#'
#' @return Returns a list formatted data.frames.
#'
#' @import dplyr
#' @import purrr
#' @noRd
get_bmd_list <- function(Enrichment_data_BMD){
  #expi_vec <- Enrichment_data_BMD |> dplyr::pull(expi) |> unique()
  expi_vec <- Enrichment_data_BMD |> dplyr::pull(Experiment) |> unique()

  # BMD enrichment: aggregate by KE (min padj)
  df_bmd_list <- expi_vec |> purrr::map(\(x){
    #Enrichment_data_BMD |> dplyr::filter(expi=={{x}}) |> 
    Enrichment_data_BMD |> dplyr::filter(Experiment=={{x}}) |> 
      dplyr::select(ke=Ke_description, BMD, BMDL, BMDU, BMD_padj=padj) |> 
      dplyr::distinct() |> 
      dplyr::group_by(ke) |> 
      dplyr::reframe(BMD=BMD |> unique() |> stringr::str_c(collapse=", "), 
        BMDL=BMDL |> unique() |> stringr::str_c(collapse=", "), 
        BMDU=BMDU |> unique() |> stringr::str_c(collapse=", "), 
        BMD_padj=BMD_padj |> unique() |> stringr::str_c(collapse=", ")
      )
  }) |> setNames(expi_vec)
  
  return(df_bmd_list=df_bmd_list)
}

#' get_docking_ke_mapping
#'
#' @description Get Docking Data
#'
#' @return Returns data.frame
#'
#' @import dplyr
#' @importFrom readxl read_excel
#' @noRd
get_docking_ke_mapping <- function(){
  # get formatted df for docking results
  MIE_gene <- readxl::read_excel("data/MIE_gene.xlsx")
  docking_matrix <- readxl::read_excel("data/docking_matrix.xlsx")
  docking_matrix <- t(docking_matrix)
  colnames(docking_matrix) <- docking_matrix[1,]
  docking_matrix <- docking_matrix[-1,]
  docking_matrix <- cbind(rownames(docking_matrix),docking_matrix)
  colnames(docking_matrix)[1] <- "GENE"
  docking_matrix <- as.data.frame(docking_matrix)
  docking_matrix <- docking_matrix |> dplyr::mutate(PFOA=PFOA |> as.numeric(), PFOS=PFOS |> as.numeric())
  KE_annotated_PFOS <- merge(docking_matrix[,c(1,3)], MIE_gene, by.x = "GENE", by.y = "GENE")
  KE_annotated_PFOA <- merge(docking_matrix[,c(1,2)], MIE_gene, by.x = "GENE", by.y = "GENE")
  colnames(KE_annotated_PFOS)[2] = colnames(KE_annotated_PFOA)[2] = "Binding"
  KE_annotated_PFOS <- cbind("PFAS" = "PFOS", KE_annotated_PFOS)
  KE_annotated_PFOA <- cbind("PFAS" = "PFOA", KE_annotated_PFOA)
  KE_annotated <- rbind(KE_annotated_PFOA, KE_annotated_PFOS)
  return(KE_annotated)
}

#' get_docking_list 
#'
#' @description Get formatted list of data.frames for each experiment from Docking Data
#'
#' @return Returns a list formatted data.frames.
#'
#' @import dplyr
#' @import purrr
#' @noRd
get_docking_list <- function(KE_annotated, expi_vec){
  #expi_vec <- Enrichment_data_BMD |> dplyr::pull(expi) |> unique()

  # Docking annotations filtered by PFAS token from experiment name
  df_docking_list <- expi_vec |> purrr::map(\(x){
    pfas_token <- strsplit(x, "_")[[1]][1];
    KE_annotated |> dplyr::filter(stringr::str_detect(PFAS, pattern=pfas_token)) |> 
      dplyr::select(PFAS, ke=Ke_description, GENE, Binding) |> 
      dplyr::distinct() |>
      dplyr::group_by(ke) |> 
      dplyr::reframe(PFAS=PFAS |> na.omit() |> unique() |> stringr::str_c(collapse=", "), 
        GENE=GENE |> na.omit() |> unique() |> stringr::str_c(collapse=", "), 
        Binding=Binding |> na.omit() |> unique() |> stringr::str_c(collapse=", ")
      )
  }) |> setNames(expi_vec)
  
  return(df_docking_list=df_docking_list)
}

#' get_paths_to_plot 
#'
#' @description Get paths to target Adverse Outcome of interest
#'
#' @return Returns updated list of nodes and edges with paths information.
#'
#' @import tibble
#' @import dplyr
#' @import purrr
#' @import igraph
#' @import tidyr
#' @import igraph
#' @export
#' @noRd
get_paths_to_plot <- function(inferred_igraph, inferred_visNet_data, ke_list_exp, ke_all_deg_list, docking_list, Enrichment_data_BMD, expi, target_node="Increased, Liver Steatosis", normalize=NULL,  PBK_time = 40){
  inferred_igraph_desc <- inferred_igraph |> igraph::set_vertex_attr(name="name", index=seq(inferred_igraph |> igraph::vcount()), value=inferred_igraph |> igraph::vertex_attr(name="ke_description", index=seq(inferred_igraph |> igraph::vcount())))
  
  if(target_node %in% V(inferred_igraph_desc)$name){
    # g_main = inferred_igraph_desc
    paths_res <- compute_paths_for_experiment(expi, inferred_igraph_desc, ke_list_exp, docking_list, Enrichment_data_BMD, ke_all_deg_list, target_node=target_node, normalize=normalize, PBK_time = PBK_time)
    
    paths_res.scores <- paths_res |> purrr::pluck("path_list") |> purrr::map("path_scores") |> tibble::enframe() |> tidyr::unnest(cols="value")
    paths_res.paths <- paths_res |> purrr::pluck("path_list") |> purrr::map("paths") |> purrr::map(\(x){x |> purrr::map(\(y){y |> names() |> stringr::str_c(collapse="; ")}) |> unlist()}) |> tibble::enframe() |> tidyr::unnest(col="value")
    
    paths_res.paths.ranked <- paths_res.paths |> dplyr::rename(path=value) |> dplyr::bind_cols(paths_res.scores |> dplyr::select(score=value)) |> dplyr::arrange(dplyr::desc(score)) |> tibble::rowid_to_column("rank") 
    paths_res.paths.ranked.nodes <- paths_res.paths.ranked |> dplyr::select(rank, path) |> tidyr::separate_longer_delim(cols="path", delim="; ") |> dplyr::rename(ke=path)
    
    inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::left_join(paths_res.paths.ranked.nodes |> dplyr::group_by(ke) |> dplyr::reframe(rank=stringr::str_c(rank, collapse=", ")), by=c("ke_description"="ke"))
    return(list("inferred_visNet_data"=inferred_visNet_data,"paths_res"=paths_res))
  }else{
    print("Target node is not present in the graph")
    return(NULL)
  }
  
}

#' Build network/KE lists per experiment
#'
#' @param probability_matrices Named list of 3D arrays (KEs x samples x timepoints) per experiment.
#' @param Enrichment_data_DEG Data frame with columns: expi, Ke_description, padj.
#' @param Enrichment_data_BMD Data frame with columns: expi, Ke_description, BMD, padj.
#' @param KE_annotated Data frame with columns: PFAS, Ke, Binding.
#' @param Biological_system_annotations Data frame with column: key_event_name (and other metadata used by make_visNetwork).
#' @param experiments Optional character vector of experiment IDs; default: names(probability_matrices).
#' @param PBK_time Integer index of the timepoint to use for the "prob" column. Default: 10.
#' @param max_path_length Integer; default 2.
#' @param n_AOs Integer; default 1.
#' @param n_MIEs Integer; default 1.
#' @param enlarge_ke_selection Logical; default TRUE.
#' @param ke_id Character; default "ke".
#' @param numerical_variables Character vector; default c("prob").
#' @param pval_variable Character; default "padj".
#' @param gene_variable Character; default "Genes".
#' @param convert_to_gene_symbols Logical; default FALSE.
#' @param group_by Character; group variable for plotting; default "aop".
#'
#' @return A named list with elements:
#' \itemize{
#'   \item net_list: per-experiment nodes/edges from make_visNetwork
#'   \item ke_list_exp: per-experiment list of KE time-series vectors
#'   \item visnet_list: per-experiment visNetwork objects
#' }
#' @export
build_network_lists <- function(
    probability_matrices,
    Enrichment_data_DEG,
    Enrichment_data_BMD,
    KE_annotated,
    experiments = names(probability_matrices),
    PBK_time = 10,
    max_path_length = 2,
    n_AOs = 1,
    n_MIEs = 1,
    mode = "out",
    enlarge_ke_selection = TRUE,
    ke_id = "ke",
    numerical_variables = c("prob"),
    pval_variable = "padj",
    gene_variable = "Genes",
    convert_to_gene_symbols = FALSE,
    group_by = "aop",
    target_node = NULL,
    normalize = NULL
) {
  
  Biological_system_annotations = AOPfingerprintR::Biological_system_annotations
  
  # Basic input checks
  if (is.null(experiments) || length(experiments) == 0) {
    stop("No experiment names found. Provide 'experiments' or ensure names(probability_matrices) are set.")
  }
  req_cols_deg <- c("Experiment", "Ke_description", "padj")
  if (!all(req_cols_deg %in% colnames(Enrichment_data_DEG))) {
    stop("Enrichment_data_DEG must contain columns: ", paste(req_cols_deg, collapse = ", "))
  }
  req_cols_ke_annot <- c("PFAS", "Ke", "Binding")
  if (!all(req_cols_ke_annot %in% colnames(KE_annotated)) ) {
    stop("KE_annotated must contain columns: PFAS, Ke, Binding")
  }
  if (!("key_event_name" %in% colnames(Biological_system_annotations))) {
    stop("Biological_system_annotations must contain 'key_event_name'.")
  }
  
  # Initialize outputs
  net_list      <- list()
  ke_list_exp   <- list()
  visnet_list   <- list()
  igraph_list   <- list()
  pal_node_fill_list <- list()
  # Optional internal list (not returned): KE present only in DEG
  ke_only_deg_list <- list()
  ke_all_deg_list <- list()
  
  docking_list = list()
  detailed_results_list <- list()
  
  optimal_path_analysis_res = list()
  
  for (expi in experiments) {
    message("Processing experiment: ", expi)
    
    # Extract 3D matrix: KEs x samples x timepoints
    mi3d <- probability_matrices[[expi]]
   
    if (is.null(mi3d)) {
      warning("No matrix found for experiment '", expi, "'. Skipping.")
      next
    }
    if (length(dim(mi3d)) != 3) {
      stop("Matrix for '", expi, "' must be 3-dimensional (KEs x samples x timepoints).")
    }
    
    n_time <- dim(mi3d)[3]
    # Build mean-over-samples for each timepoint (KEs x timepoints)
    mean_over_samples <- NULL
    for (tp in seq_len(n_time)) {
      # Row means across samples for timepoint tp
      rms <- rowMeans(mi3d[, , tp], na.rm = TRUE)
      mean_over_samples <- cbind(mean_over_samples, rms)
    }
    mi <- mean_over_samples
    rownames(mi) <- rownames(mi3d) # ensure KE names present
    
    # Per-KE time series list
    ke_list <- vector(mode = "list", length = nrow(mi))
    names(ke_list) <- rownames(mi)
    for (i in seq_len(nrow(mi))) {
      ke_list[[rownames(mi)[i]]] <- mi[i, ]
    }
    
    ke_list_exp[[expi]]  <- ke_list
    
    # Base DF with prob at PBK_time
    if (PBK_time < 1 || PBK_time > ncol(mi)) {
      stop("PBK_time=", PBK_time, " is out of range. Available timepoints: 1..", ncol(mi))
    }
    df <- data.frame(
      ke   = rownames(mi),
      PBK_prob = mi[, PBK_time],
      stringsAsFactors = FALSE
    )
    # DEG enrichment: aggregate by KE (min padj)
    df_deg <- Enrichment_data_DEG[Enrichment_data_DEG$Experiment == expi, c("Ke_description","relevantGenesInGeneSet", "padj","log2FoldChange")]
    
    df_deg <- dplyr::as_tibble(df_deg) |>
      dplyr::group_by(Ke_description) |>
      dplyr::slice_min(order_by=padj, n=1, with_ties=FALSE) |> dplyr::ungroup() |>
      as.data.frame()
    
    df_deg$pi_value = -log10(as.numeric(df_deg$padj)) * abs(as.numeric(df_deg$log2FoldChange))
    df_deg <- df_deg |> dplyr::mutate(prob=pi_value)
    
    ke_all_deg_list[[expi]] = df_deg
    # NOTE: reuse 'prob' column for padj to align with downstream merge
    df_deg <- df_deg |> dplyr::select(ke=Ke_description, pi_value, DE_padj=padj, log2FoldChange) |> dplyr::mutate(log2FoldChange=log2FoldChange |> as.numeric()) |> dplyr::distinct()
  
    # BMD enrichment: aggregate by KE (min padj)
    df_bmd <- Enrichment_data_BMD |> dplyr::filter(Experiment=={{expi}}) |> dplyr::select(ke=Ke_description, BMD, BMDL, BMDU, BMD_padj=padj) |> dplyr::distinct()
    df_bmd <- df_bmd |> dplyr::group_by(ke) |> dplyr::reframe(BMD=BMD |> unique() |> stringr::str_c(collapse=", "), BMDL=BMDL |> unique() |> stringr::str_c(collapse=", "), BMDU=BMDU |> unique() |> stringr::str_c(collapse=", "), BMD_padj=BMD_padj |> unique() |> stringr::str_c(collapse=", "))

    # Docking annotations filtered by PFAS token from experiment name
    pfas_token <- strsplit(expi, "_")[[1]][1]
    df_docking <- KE_annotated |> dplyr::filter(stringr::str_detect(PFAS, pattern=pfas_token)) |> dplyr::select(PFAS, ke=Ke_description, GENE, Binding) |> dplyr::distinct()
    docking_list[[expi]] = df_docking
    df_docking <- df_docking |> dplyr::group_by(ke) |> dplyr::reframe(PFAS=PFAS |> na.omit() |> unique() |> stringr::str_c(collapse=", "), GENE=GENE |> na.omit() |> unique() |> stringr::str_c(collapse=", "), Binding=Binding |> na.omit() |> unique() |> stringr::str_c(collapse=", "))
    
    detailed_results <- df |> dplyr::full_join(df_deg, by=c("ke"="ke")) |> dplyr::full_join(df_bmd, by=c("ke"="ke")) |> dplyr::full_join(df_docking, by=c("ke"="ke")) |> dplyr::left_join(AOPfingerprintR::Biological_system_annotations, by=c("ke"="key_event_name"))
    
    # Add required fields
    detailed_results <- detailed_results |> dplyr::mutate(Experiment=expi) |> dplyr::rename(Ke_description=ke, ke=ke.y)
    
    ke_desc_to_ke_id_map = unique(AOPfingerprintR::aop_ke_table_hure[,c("Ke","Ke_description")])
    rownames(ke_desc_to_ke_id_map) = ke_desc_to_ke_id_map$Ke_description
    
    idx = which(is.na(detailed_results$ke))
    if(length(idx)>0){
      detailed_results$ke[idx] =  ke_desc_to_ke_id_map[detailed_results$Ke_description[idx],"Ke"]
      
    }
    
    # Get igraph
    ke_set <- detailed_results |> dplyr::pull({{ke_id}}) |> unique()
    inferred_igraph <- get_igraph(ke_set=ke_set, 
                                  enlarge_ke_selection=enlarge_ke_selection, 
                                  max_path_length=max_path_length, 
                                  n_AOs=n_AOs, 
                                  n_MIEs=n_MIEs, 
                                  mode=mode)
    inferred_visNet_data <- visNetwork::toVisNetworkData(inferred_igraph)

    # address multiple mapping of KE, AO, and MIE
    aop_ke_df <- AOPfingerprintR::aop_ke_table_hure |> dplyr::select(Ke, Ke_type) |> dplyr::add_count(Ke_type) |> dplyr::distinct() |> tidyr::pivot_wider(names_from=Ke_type, values_from=n) |> dplyr::mutate(Ke_type=ifelse(is.na(AdverseOutcome),ifelse(is.na(MolecularInitiatingEvent),"KeyEvent","MolecularInitiatingEvent"),ifelse(is.na(MolecularInitiatingEvent),"AdverseOutcome",ifelse(AdverseOutcome>MolecularInitiatingEvent,"AdverseOutcome","MolecularInitiatingEvent")))) 

    # add node shape by KE type
    inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::left_join(aop_ke_df |> dplyr::select(Ke, ke_type=Ke_type), by=c("id"="Ke")) |> dplyr::relocate(ke_type, .after=ke_description)
    node_shapes <- setNames(object=c("triangle","square","star"), nm=c("MolecularInitiatingEvent","KeyEvent","AdverseOutcome"))
    inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::mutate(shape=node_shapes[ke_type] |> unname())

    # format organ info for AOP
    Annotate_AOPs <- AOPfingerprintR::Annotate_AOPs |> dplyr::mutate(Organ=ifelse(is.na(Organ),"Not annotated to organ",Organ), SSbD_category_organ=paste0(SSbD_category, " - ", Organ)) |> dplyr::distinct()

    # add vis groups
    #inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::left_join(AOPfingerprintR::aop_ke_table_hure |> dplyr::select(Ke, Aop, aop_name=a.name), by=c("id"="Ke")) |> dplyr::distinct() |> dplyr::left_join(Annotate_AOPs |> dplyr::select(AOP, SSbD_category_organ), by=c("Aop"="AOP")) |> dplyr::select(-Aop)  |> dplyr::group_by(id) |> dplyr::summarize(id, ke_description, ke_type, label, shape, aop=paste0(aop_name |> na.omit() |> unique() |> stringr::str_c(collapse=", ")), ssbd=paste0(SSbD_category_organ |> na.omit() |> unique() |> stringr::str_c(collapse=", "))) |> dplyr::ungroup() |> dplyr::distinct()
    inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::left_join(AOPfingerprintR::aop_ke_table_hure |> dplyr::select(Ke, Aop, aop_name=a.name), by=c("id"="Ke")) |> dplyr::distinct() |> dplyr::left_join(Annotate_AOPs |> dplyr::select(AOP, SSbD_category_organ), by=c("Aop"="AOP")) |> dplyr::select(-Aop) |> dplyr::group_by(id) |> dplyr::reframe(id, ke_description, ke_type, label, shape, aop=aop_name |> na.omit() |> unique() |> paste0(collapse=", "), ssbd=SSbD_category_organ |> na.omit() |> unique() |> paste0(collapse=", ")) |> dplyr::ungroup() |> dplyr::distinct()

    # get paths to target AO
    if(!is.null(target_node)){
      path_analysis_res <- get_paths_to_plot(inferred_igraph, inferred_visNet_data, ke_list_exp, ke_all_deg_list, docking_list, Enrichment_data_BMD, expi, target_node=target_node, normalize=normalize,PBK_time=PBK_time)
      if(!is.null(path_analysis_res)){
        inferred_visNet_data = path_analysis_res$inferred_visNet_data
        optimal_path_analysis_res[[expi]] = path_analysis_res$paths_res
      }
    }

    # add all scores
    inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::left_join(detailed_results |> dplyr::select(ke, PBK_prob, pi_value, DE_padj, log2FoldChange, BMD, BMDL, BMDU, BMD_padj, PFAS, GENE, Binding), by=c("id"="ke")) |> dplyr::distinct()

    # fix label names
    inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::mutate(label=ke_description)
    inferred_visNet_data$edges <- inferred_visNet_data$edges |> dplyr::mutate(label=type.r)

    # add vis display text
    inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::mutate(title=paste0("<p> Id:", id, "</p>", "<p> Description:", label, "</p>"))
    inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::mutate(title=paste0(title,"<p> AOPs:", aop, "</p>", "<p> Hazard Class:", ssbd, "</p>"))
    score_cols <- c("PBK_prob", "pi_value", "DE_padj", "log2FoldChange", "BMD", "BMDL", "BMDU", "BMD_padj", "PFAS", "GENE", "Binding")
    for(x in score_cols){
      inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::mutate(title=paste0(title, "<p> ",x,":", get(x), "</p>"))
    }

    # evidence counts
    inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::rowwise() |> dplyr::mutate(evidence=c(PBK_prob, DE_padj, Binding) |> is.na() |> magrittr::not() |> sum())

    # node fill color 
    my_blues <- colorRampPalette(c(
      "skyblue",
      "navyblue" # custom starting blue for low values
    ))

    na.color <- "white" 
    pal_node_fill <- scales::col_numeric(
      palette = my_blues(100),
      domain  = range(c(0,inferred_visNet_data$nodes$PBK_prob,100), na.rm = TRUE),
      na.color = na.color
    )
    
    # inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::mutate(color.background=PBK_prob |> pal_node_fill())
    
    zero_prob_color <- "#BFBFBF"
    
    inferred_visNet_data$nodes <- inferred_visNet_data$nodes |>
      dplyr::mutate(
        color.background = dplyr::case_when(
          is.na(PBK_prob) ~ na.color,
          PBK_prob == 0 ~ zero_prob_color,
          TRUE ~ pal_node_fill(PBK_prob)
        )
      )
    
    # node border color
    pal_node_border <- c("black", "seagreen", "#D4A017", "red") |> setNames(0:3)
    inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::mutate(color.border=pal_node_border[evidence |> as.character()] |> unname())

    # node border width
    inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::mutate(borderWidth=ifelse(evidence==0, 0.5, 5))

    # remove edge labels
    inferred_visNet_data$edges = inferred_visNet_data$edges[,1:2]
    inferred_visNet_data$edges = unique(inferred_visNet_data$edges)

    # Get visNetwork
    # nodes=inferred_visNet_data$nodes
    # edges=inferred_visNet_data$edges
    # group_by=group_by
    # col_pal = pal_node_fill
    vn <- get_visnet(nodes=inferred_visNet_data$nodes, edges=inferred_visNet_data$edges, group_by=group_by, col_pal = pal_node_fill)

    # Collect outputs
    igraph_list[[expi]] <- inferred_igraph
    visnet_list[[expi]]  <- vn
    net_list[[expi]]     <- inferred_visNet_data
    pal_node_fill_list[[expi]]     <- pal_node_fill
    detailed_results_list[[expi]] <- detailed_results
  }
  
  # Return the three lists as requested
  list(
    net_list     = net_list,
    ke_list_exp  = ke_list_exp,
    visnet_list  = visnet_list,
    igraph_list  = igraph_list,
    pal_node_fill_list  = pal_node_fill_list,
    ke_only_deg_list = ke_only_deg_list,
    docking_list = docking_list,
    ke_all_deg_list = ke_all_deg_list,
    detailed_results_list = detailed_results_list,
    optimal_path_analysis_res = optimal_path_analysis_res
  )
}


#' Build network/KE lists per experiment with only PBK BMD
#'
#' @param probability_matrices Named list of 3D arrays (KEs x samples x timepoints) per experiment.
#' @param Enrichment_data_DEG Data frame with columns: expi, Ke_description, padj.
#' @param Enrichment_data_BMD Data frame with columns: expi, Ke_description, BMD, padj.
#' @param KE_annotated Data frame with columns: PFAS, Ke, Binding.
#' @param Biological_system_annotations Data frame with column: key_event_name (and other metadata used by make_visNetwork).
#' @param experiments Optional character vector of experiment IDs; default: names(probability_matrices).
#' @param PBK_time Integer index of the timepoint to use for the "prob" column. Default: 10.
#' @param max_path_length Integer; default 2.
#' @param n_AOs Integer; default 1.
#' @param n_MIEs Integer; default 1.
#' @param enlarge_ke_selection Logical; default TRUE.
#' @param ke_id Character; default "ke".
#' @param numerical_variables Character vector; default c("prob").
#' @param pval_variable Character; default "padj".
#' @param gene_variable Character; default "Genes".
#' @param convert_to_gene_symbols Logical; default FALSE.
#' @param group_by Character; group variable for plotting; default "aop".
#'
#' @return A named list with elements:
#' \itemize{
#'   \item net_list: per-experiment nodes/edges from make_visNetwork
#'   \item ke_list_exp: per-experiment list of KE time-series vectors
#'   \item visnet_list: per-experiment visNetwork objects
#' }
#' @export
build_network_lists_only_PBK_BMD <- function(
    probability_matrices,
    Enrichment_data_BMD,
    experiments = names(probability_matrices),
    PBK_time = 10,
    max_path_length = 2,
    n_AOs = 1,
    n_MIEs = 1,
    mode = "out",
    enlarge_ke_selection = TRUE,
    ke_id = "ke",
    numerical_variables = c("prob"),
    pval_variable = "padj",
    gene_variable = "Genes",
    convert_to_gene_symbols = FALSE,
    group_by = "aop",
    normalize = NULL
) {
  Biological_system_annotations = AOPfingerprintR::Biological_system_annotations
  # Basic input checks
  if (is.null(experiments) || length(experiments) == 0) {
    stop("No experiment names found. Provide 'experiments' or ensure names(probability_matrices) are set.")
  }
  
  if (!("key_event_name" %in% colnames(Biological_system_annotations))) {
    stop("Biological_system_annotations must contain 'key_event_name'.")
  }
  
  # Initialize outputs
  net_list      <- list()
  ke_list_exp   <- list()
  visnet_list   <- list()
  
  for (expi in experiments) {
    message("Processing experiment: ", expi)
    
    # Extract 3D matrix: KEs x samples x timepoints
    mi3d <- probability_matrices[[expi]]
    if (is.null(mi3d)) {
      warning("No matrix found for experiment '", expi, "'. Skipping.")
      next
    }
    if (length(dim(mi3d)) != 3) {
      stop("Matrix for '", expi, "' must be 3-dimensional (KEs x samples x timepoints).")
    }
    
    n_time <- dim(mi3d)[3]
    # Build mean-over-samples for each timepoint (KEs x timepoints)
    mean_over_samples <- NULL
    for (tp in seq_len(n_time)) {
      # Row means across samples for timepoint tp
      rms <- rowMeans(mi3d[, , tp], na.rm = TRUE)
      mean_over_samples <- cbind(mean_over_samples, rms)
    }
    mi <- mean_over_samples
    rownames(mi) <- rownames(mi3d) # ensure KE names present
    
    # Per-KE time series list
    ke_list <- vector(mode = "list", length = nrow(mi))
    names(ke_list) <- rownames(mi)
    for (i in seq_len(nrow(mi))) {
      ke_list[[rownames(mi)[i]]] <- mi[i, ]
    }
    
    # Base DF with prob at PBK_time
    if (PBK_time < 1 || PBK_time > ncol(mi)) {
      stop("PBK_time=", PBK_time, " is out of range. Available timepoints: 1..", ncol(mi))
    }
    df <- data.frame(
      ke   = rownames(mi),
      PBK_prob = mi[, PBK_time],
      stringsAsFactors = FALSE
    )
    
    # BMD enrichment: aggregate by KE (min padj)
    df_bmd <- Enrichment_data_BMD |> dplyr::filter(Experiment=={{expi}}) |> dplyr::select(ke=Ke_description, BMD, BMDL, BMDU, BMD_padj=padj) |> dplyr::distinct()
    df_bmd <- df_bmd |> dplyr::group_by(ke) |> dplyr::reframe(BMD=BMD |> unique() |> stringr::str_c(collapse=", "), BMDL=BMDL |> unique() |> stringr::str_c(collapse=", "), BMDU=BMDU |> unique() |> stringr::str_c(collapse=", "), BMD_padj=BMD_padj |> unique() |> stringr::str_c(collapse=", "))
    
    # Docking annotations filtered by PFAS token from experiment name
    pfas_token <- strsplit(expi, "_")[[1]][1]
    detailed_results <- df |> dplyr::full_join(df_bmd, by=c("ke"="ke")) |> dplyr::left_join(AOPfingerprintR::Biological_system_annotations, by=c("ke"="key_event_name"))
    
    # Add required fields
    detailed_results <- detailed_results |> dplyr::mutate(Experiment=expi) |> dplyr::rename(Ke_description=ke, ke=ke.y)
    
    # Get igraph
    ke_set <- detailed_results |> dplyr::pull({{ke_id}}) |> unique()
    inferred_igraph <- get_igraph(ke_set=ke_set, 
                                  enlarge_ke_selection=enlarge_ke_selection, 
                                  max_path_length=max_path_length, 
                                  n_AOs=n_AOs, 
                                  n_MIEs=n_MIEs, 
                                  mode=mode)
    inferred_visNet_data <- visNetwork::toVisNetworkData(inferred_igraph)
    
    # address multiple mapping of KE, AO, and MIE
    aop_ke_df <- AOPfingerprintR::aop_ke_table_hure |> dplyr::select(Ke, Ke_type) |> dplyr::add_count(Ke_type) |> dplyr::distinct() |> tidyr::pivot_wider(names_from=Ke_type, values_from=n) |> dplyr::mutate(Ke_type=ifelse(is.na(AdverseOutcome),ifelse(is.na(MolecularInitiatingEvent),"KeyEvent","MolecularInitiatingEvent"),ifelse(is.na(MolecularInitiatingEvent),"AdverseOutcome",ifelse(AdverseOutcome>MolecularInitiatingEvent,"AdverseOutcome","MolecularInitiatingEvent")))) 
    
    # add node shape by KE type
    inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::left_join(aop_ke_df |> dplyr::select(Ke, ke_type=Ke_type), by=c("id"="Ke")) |> dplyr::relocate(ke_type, .after=ke_description)
    node_shapes <- setNames(object=c("triangle","square","star"), nm=c("MolecularInitiatingEvent","KeyEvent","AdverseOutcome"))
    inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::mutate(shape=node_shapes[ke_type] |> unname())
    
    # format organ info for AOP
    Annotate_AOPs <- AOPfingerprintR::Annotate_AOPs |> dplyr::mutate(Organ=ifelse(is.na(Organ),"Not annotated to organ",Organ), SSbD_category_organ=paste0(SSbD_category, " - ", Organ)) |> dplyr::distinct()
    
    # add vis groups
    #inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::left_join(AOPfingerprintR::aop_ke_table_hure |> dplyr::select(Ke, Aop, aop_name=a.name), by=c("id"="Ke")) |> dplyr::distinct() |> dplyr::left_join(Annotate_AOPs |> dplyr::select(AOP, SSbD_category_organ), by=c("Aop"="AOP")) |> dplyr::select(-Aop)  |> dplyr::group_by(id) |> dplyr::summarize(id, ke_description, ke_type, label, shape, aop=paste0(aop_name |> na.omit() |> unique() |> stringr::str_c(collapse=", ")), ssbd=paste0(SSbD_category_organ |> na.omit() |> unique() |> stringr::str_c(collapse=", "))) |> dplyr::ungroup() |> dplyr::distinct()
    inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::left_join(AOPfingerprintR::aop_ke_table_hure |> dplyr::select(Ke, Aop, aop_name=a.name), by=c("id"="Ke")) |> dplyr::distinct() |> dplyr::left_join(Annotate_AOPs |> dplyr::select(AOP, SSbD_category_organ), by=c("Aop"="AOP")) |> dplyr::select(-Aop) |> dplyr::group_by(id) |> dplyr::reframe(id, ke_description, ke_type, label, shape, aop=aop_name |> na.omit() |> unique() |> paste0(collapse=", "), ssbd=SSbD_category_organ |> na.omit() |> unique() |> paste0(collapse=", ")) |> dplyr::ungroup() |> dplyr::distinct()
    
    # get paths to target AO
    ### write a wrapper for this
    inferred_igraph_desc <- inferred_igraph |> igraph::set_vertex_attr(name="name", index=seq(inferred_igraph |> igraph::vcount()), value=inferred_igraph |> igraph::vertex_attr(name="ke_description", index=seq(inferred_igraph |> igraph::vcount())))
    
    # add all scores
    inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::left_join(detailed_results |> dplyr::select(ke, PBK_prob, BMD, BMDL, BMDU), by=c("id"="ke")) |> dplyr::distinct()
    
    # fix label names
    inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::mutate(label=ke_description)
    inferred_visNet_data$edges <- inferred_visNet_data$edges |> dplyr::mutate(label=type.r)
    
    # add vis display text
    inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::mutate(title=paste0("<p> Id:", id, "</p>", "<p> Description:", label, "</p>"))
    inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::mutate(title=paste0(title,"<p> AOPs:", aop, "</p>", "<p> Hazard Class:", ssbd, "</p>"))
    score_cols <- c("PBK_prob", "BMD", "BMDL", "BMDU")
    for(x in score_cols){
      inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::mutate(title=paste0(title, "<p> ",x,":", get(x), "</p>"))
    }
    
    #inferred_visNet_data$nodes$evidence = inferred_visNet_data$nodes$PBK_prob
    #inferred_visNet_data$nodes$evidence[is.na(inferred_visNet_data$nodes$evidence)] = 0
    #
    #library(dplyr)
    #library(scales)
    
    # node fill color 
    my_blues <- colorRampPalette(c(
      "skyblue",
      "navyblue" # custom starting blue for low values
    ))

    na.color <- "white" 
    pal_node_fill <- scales::col_numeric(
      palette = my_blues(100),
      domain  = range(c(0,inferred_visNet_data$nodes$PBK_prob,100), na.rm = TRUE),
      na.color = na.color
    )
    
    inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::mutate(color.background=PBK_prob |> pal_node_fill())

    ## Continuous palette for positive values
    #pos_col_fun <- col_numeric(
    #  palette = c("seagreen", "orange", "red"),
    #  domain  = range(inferred_visNet_data$nodes$evidence[
    #    inferred_visNet_data$nodes$evidence > 0
    #  ],
    #  na.rm = TRUE)
    #)
    #
    #inferred_visNet_data$nodes <- inferred_visNet_data$nodes %>%
    #  mutate(color = ifelse(
    #    evidence == 0,
    #    "lightblue",           # fixed color for zero
    #    pos_col_fun(evidence)  # continuous color for >0
    #  ))
    
    # node border color
    inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::mutate(color.border="black")

    # node border width
    inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::mutate(borderWidth=1)

    # remove edge labels
    inferred_visNet_data$edges = inferred_visNet_data$edges[,1:2]
    inferred_visNet_data$edges = unique(inferred_visNet_data$edges)

    # Get visNetwork
    vn <- get_visnet(nodes=inferred_visNet_data$nodes, edges=inferred_visNet_data$edges, group_by=group_by, col_pal=pal_node_fill)
    
    # Collect outputs
    visnet_list[[expi]]  <- vn
    net_list[[expi]]     <- inferred_visNet_data
    ke_list_exp[[expi]]  <- ke_list
  }
  
  # Return the three lists as requested
  list(
    net_list     = net_list,
    ke_list_exp  = ke_list_exp,
    visnet_list  = visnet_list
  )
}

#' Build network/KE lists per experiment
#'
#' @param reverse dosimetry summarized results per experiment.
#' @param Biological_system_annotations Data frame with column: key_event_name (and other metadata used by make_visNetwork).
#' @param experiments Optional character vector of experiment IDs; default: names(probability_matrices).
#' @param max_path_length Integer; default 2.
#' @param n_AOs Integer; default 1.
#' @param n_MIEs Integer; default 1.
#' @param enlarge_ke_selection Logical; default TRUE.
#' @param ke_id Character; default "ke".
#' @param group_by Character; group variable for plotting; default "aop".
#'
#' @return A named list with elements:
#' \itemize{
#'   \item net: per-experiment nodes/edges from make_visNetwork
#'   \item visnet: per-experiment visNetwork objects
#' }
#'
#' @import dplyr
#' @import tidyr
#' @import scales
#' @import AOPfingerprintR
#'
#' @export
get_ke_net_reverse_dosimetry <- function(
    reverse_dosimetry_summarized,
    max_path_length = 2,
    n_AOs = 1,
    n_MIEs = 1,
    mode = "out",
    enlarge_ke_selection = TRUE,
    ke_id = "ke",
    group_by = "aop",
    normalize = NULL
){
  Biological_system_annotations <- AOPfingerprintR::Biological_system_annotations
  DF <- reverse_dosimetry_summarized |> dplyr::mutate(exposure=exposure |> round())
  DF <- DF |> dplyr::inner_join(Biological_system_annotations |> dplyr::select(ke, key_event_name) |> dplyr::distinct(), by=c("key_event"="key_event_name"))

  # Get igraph
  ke_set <- DF |> dplyr::pull({{ke_id}}) |> unique()
  inferred_igraph <- get_igraph(ke_set=ke_set, 
                                enlarge_ke_selection=enlarge_ke_selection, 
                                max_path_length=max_path_length, 
                                n_AOs=n_AOs, 
                                n_MIEs=n_MIEs, 
                                mode=mode)
  inferred_visNet_data <- visNetwork::toVisNetworkData(inferred_igraph)
  
  # address multiple mapping of KE, AO, and MIE
  aop_ke_df <- AOPfingerprintR::aop_ke_table_hure |> dplyr::select(Ke, Ke_type) |> dplyr::add_count(Ke_type) |> dplyr::distinct() |> tidyr::pivot_wider(names_from=Ke_type, values_from=n) |> dplyr::mutate(Ke_type=ifelse(is.na(AdverseOutcome),ifelse(is.na(MolecularInitiatingEvent),"KeyEvent","MolecularInitiatingEvent"),ifelse(is.na(MolecularInitiatingEvent),"AdverseOutcome",ifelse(AdverseOutcome>MolecularInitiatingEvent,"AdverseOutcome","MolecularInitiatingEvent")))) 

  # add node shape by KE type
  inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::left_join(aop_ke_df |> dplyr::select(Ke, ke_type=Ke_type), by=c("id"="Ke")) |> dplyr::relocate(ke_type, .after=ke_description)
  node_shapes <- setNames(object=c("triangle","square","star"), nm=c("MolecularInitiatingEvent","KeyEvent","AdverseOutcome"))
  inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::mutate(shape=node_shapes[ke_type] |> unname())

  # format organ info for AOP
  Annotate_AOPs <- AOPfingerprintR::Annotate_AOPs |> dplyr::mutate(Organ=ifelse(is.na(Organ),"Not annotated to organ",Organ), SSbD_category_organ=paste0(SSbD_category, " - ", Organ)) |> dplyr::distinct()

  # add vis groups
  #inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::left_join(AOPfingerprintR::aop_ke_table_hure |> dplyr::select(Ke, Aop, aop_name=a.name), by=c("id"="Ke")) |> dplyr::distinct() |> dplyr::left_join(Annotate_AOPs |> dplyr::select(AOP, SSbD_category_organ), by=c("Aop"="AOP")) |> dplyr::select(-Aop)  |> dplyr::group_by(id) |> dplyr::summarize(id, ke_description, ke_type, label, shape, aop=paste0(aop_name |> na.omit() |> unique() |> stringr::str_c(collapse=", ")), ssbd=paste0(SSbD_category_organ |> na.omit() |> unique() |> stringr::str_c(collapse=", "))) |> dplyr::ungroup() |> dplyr::distinct()
  inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::left_join(AOPfingerprintR::aop_ke_table_hure |> dplyr::select(Ke, Aop, aop_name=a.name), by=c("id"="Ke")) |> dplyr::distinct() |> dplyr::left_join(Annotate_AOPs |> dplyr::select(AOP, SSbD_category_organ), by=c("Aop"="AOP")) |> dplyr::select(-Aop) |> dplyr::group_by(id) |> dplyr::reframe(id, ke_description, ke_type, label, shape, aop=aop_name |> na.omit() |> unique() |> paste0(collapse=", "), ssbd=SSbD_category_organ |> na.omit() |> unique() |> paste0(collapse=", ")) |> dplyr::ungroup() |> dplyr::distinct()

  # add all scores
  inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::left_join(DF |> dplyr::select(ke, exposure, rel_error, POD, standard_deviation=sd, variance=var), by=c("id"="ke")) |> dplyr::distinct()

  # fix label names
  inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::mutate(label=ke_description)

  # add vis display text
  inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::mutate(title=paste0("<p> Id:", id, "</p>", "<p> Description:", label, "</p>"))
  inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::mutate(title=paste0(title,"<p> AOPs:", aop, "</p>", "<p> Hazard Class:", ssbd, "</p>"))
  #score_cols <- c("PBK_prob", "pi_value", "DE_padj", "log2FoldChange", "BMD", "BMDL", "BMDU", "BMD_padj", "PFAS", "GENE", "Binding")
  score_cols <- c("exposure", "rel_error", "POD", "standard_deviation", "variance")
  for(x in score_cols){
    inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::mutate(title=paste0(title, "<p> ",x,":", get(x), "</p>"))
  }

  # node fill color 
  my_blues <- colorRampPalette(c(
    "navyblue",
    "skyblue" # custom starting blue for low values
  ))

  na.color <- "white" 
  pal_node_fill <- scales::col_numeric(
    palette = my_blues(100),
    domain  = range(inferred_visNet_data$nodes$exposure, na.rm = TRUE),
    na.color = na.color
  )
  
  inferred_visNet_data$nodes <- inferred_visNet_data$nodes |> dplyr::mutate(color.background=exposure |> pal_node_fill())

  # Get visNetwork
  vn <- get_visnet_reverse_dosimetry(nodes=inferred_visNet_data$nodes, edges=inferred_visNet_data$edges, group_by=group_by, col_pal = pal_node_fill)

  # Return the three lists as requested
  list(
    net     = inferred_visNet_data,
    visnet  = vn,
    igraph  = inferred_igraph,
    pal_node_fill  = pal_node_fill
  )
}


