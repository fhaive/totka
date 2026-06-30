
# INPUTS:
# PFAS: a string that can be either PFOS or PFOA
# tissue_of_interest: a string representing the PBK compartmnet of interest (e.g. liver, serum)
# pod_variable: a stringthat represent the point of departure value (e.g. BMD)
# time_point_invitro: an integer representing the time point of the invitro data
# tpi: integer, one of the exposure time of PBK simulation
# Simulation_results: nested list containing the PBK simulation results
# Enrichment_data: a list of KE enrichment dataframes, one for each PFAS
# enlarge_ke_to_uncategorized: if true, it includes also KE that do not have a tissue label


#' Integrate Omics KE Enrichment with PBPK Simulations for PFAS
#'
#' @description
#' Given KE (Key Event) enrichment results stratified by in vitro time points
#' and PBPK (physiologically based pharmacokinetic) model predictions for
#' internal concentrations in a target tissue/compartment, this function
#' determines which KEs are *triggered* at specified in vitro and in vivo time
#' points for one or more PFAS (PFOS, PFOA).
#' @param PFAS Character scalar. Which PFAS to process. One of:
#'   \itemize{
#'     \item \code{"PFOS"}
#'     \item \code{"PFOA"}
#'     \item \code{"All"} — process both PFOS and PFOA
#'   }
#'
#' @param tissue_of_interest Character scalar. Name of the PBPK compartment/tissue
#'   to evaluate (e.g., \code{"liver"}, \code{"serum"}). Must match a key in
#'   \code{Simulation_results[[PFAS]]}.
#'
#' @param pod_variable Character scalar. Point-of-departure variable used in
#'   the KE-triggering logic (e.g., \code{"BMD"}).
#'
#' @param time_point_invitro Integer scalar. In vitro time point used to select
#'   the corresponding KE enrichment subset (must be present in
#'   \code{Enrichment_data[[PFAS]]$time}).
#'
#' @param time_point_invivo Integer scalar. Index of the PBPK simulation time
#'   to evaluate (passed to downstream logic as \code{scaled_time_idx}).
#'
#' @param Simulation_results Nested list of PBPK predictions. Accessed as
#'   \code{Simulation_results[[PFAS]][[tissue_of_interest]]}, which must return
#'   a data frame where:
#'   \itemize{
#'     \item Column 1 is an index (e.g., simulation time).
#'     \item Columns 2+ are simulated internal concentrations for individual PBPK models
#'           and/or BMA draws (units typically \eqn{\mu g/L}).
#'   }
#'
#' @param Enrichment_data List of KE enrichment data frames, named by PFAS
#'   (e.g., \code{"PFOS"}, \code{"PFOA"}). Each data frame should include:
#'   \itemize{
#'     \item \code{organ_tissue} (character; may contain \code{NA})
#'     \item \code{time} (integer/numeric; in vitro time point)
#'     \item other KE enrichment fields used by your pipeline
#'   }
#'
#' @param enlarge_ke_to_uncategorized Logical (default \code{TRUE}). If \code{TRUE},
#'   KEs with missing \code{organ_tissue} are assigned to \code{tissue_of_interest}.
#'
#' @return
#' A list with one element:
#' \describe{
#'   \item{\code{triggered_ke_list}}{Named by PFAS. Each value is the
#'   \code{triggered_KE} list produced by
#'   \code{identify_kes_triggered_by_exposure_simulation()}, where each element
#'   is a data frame of KEs triggered for one PBPK simulation (individual model
#'   or BMA draw).}
#' }
#'
#' @export
omics_pbpk_integration = function(PFAS,
                                  tissue_of_interest,
                                  pod_variable,
                                  time_point_invitro,
                                  time_point_invivo,
                                  Simulation_results,
                                  Enrichment_data, 
                                  enlarge_ke_to_uncategorized = T){
  
  
  # ---- Input checks ----------------------------------------------------------
  #  type/length checks
  stopifnot(is.character(PFAS), length(PFAS) == 1L, !is.na(PFAS))
  stopifnot(is.character(tissue_of_interest), length(tissue_of_interest) == 1L, nzchar(tissue_of_interest))
  stopifnot(is.character(pod_variable), length(pod_variable) == 1L, nzchar(pod_variable), !is.na(pod_variable))
  stopifnot(length(time_point_invitro) == 1L, is.numeric(time_point_invitro), !is.na(time_point_invitro))
  stopifnot(length(time_point_invivo) == 1L, is.numeric(time_point_invivo), !is.na(time_point_invivo))
  stopifnot(is.list(Simulation_results))
  stopifnot(is.list(Enrichment_data))
  stopifnot(is.logical(enlarge_ke_to_uncategorized), length(enlarge_ke_to_uncategorized) == 1L, !is.na(enlarge_ke_to_uncategorized))
  
  # Allowed PFAS
  allowed_pfas <- c("PFOS","PFOA","All")
  if (!PFAS %in% allowed_pfas) {
    stop(sprintf("`PFAS` must be one of %s; got: %s",
                 paste(shQuote(allowed_pfas), collapse = ", "), PFAS))
  }
  
  # Expand PFAS = "ALL"
  PFAS_vec <- if (identical(PFAS, "ALL")) c("PFOA", "PFOS") else PFAS
  
  
  # Ensure required PFAS entries exist in both lists
  missing_enrich <- setdiff(PFAS_vec, names(Enrichment_data))
  if (length(missing_enrich)) {
    stop(sprintf("`Enrichment_data` must contain entries for: %s",
                 paste(shQuote(missing_enrich), collapse = ", ")))
  }
  missing_sim <- setdiff(PFAS_vec, names(Simulation_results))
  if (length(missing_sim)) {
    stop(sprintf("`Simulation_results` must contain entries for: %s",
                 paste(shQuote(missing_sim), collapse = ", ")))
  }
  
  # ---- Constants / metadata --------------------------------------------------
  PFOS_mw = 500.13 #g/mol 
  PFOA_mw = 414.07 #g/mol 
  PFAS_mw = c("PFOS" = PFOS_mw, "PFOA"= PFOA_mw)
  BMD_internal_dose_unit = "µM"
  PBPK_tissue_dose_unit = "µg/L"
  # unique_tp = c(1,4,10,14) # TO DO: to be taken from the input file, or to be asked to the user when uploading the KE & POD file


  # ---- Main loop -------------------------------------------------------------
  triggered_ke_list = vector("list", length(PFAS))
  

  for(ss in PFAS){
    
    # --- Validate Enrichment_data structure for this PFAS
    KE_annotated <- Enrichment_data[[ss]]
    
    if (!is.data.frame(KE_annotated)) {
      stop(sprintf("`Enrichment_data[['%s']]` must be a data.frame.", ss))
    }
    req_cols <- c("organ_tissue", "time")
    missing_cols <- setdiff(req_cols, colnames(KE_annotated))
    if (length(missing_cols)) {
      stop(sprintf("`Enrichment_data[['%s']]` is missing required columns: %s",
                   ss, paste(shQuote(missing_cols), collapse = ", ")))
    }
    if (!is.numeric(KE_annotated$time)) {
      stop(sprintf("`Enrichment_data[['%s']]$time` must be numeric/integer.", ss))
    }
    # Ensure in vitro time point is present
    if (!(time_point_invitro %in% unique(KE_annotated$time))) {
      stop(sprintf("`time_point_invitro = %s` not found in `Enrichment_data[['%s']]$time`.",
                   as.character(time_point_invitro), ss))
    }
    
    
    #Align vocabularies between AOP annotaiton and serum
    serum_related_organs = c("artery","blood","blood plasma",
                             "blood plasma/liver",
                             "blood serum",
                             "blood serum/anterior pituitary gland",
                             "blood vessel",
                             "blood vessel endothelium",
                             "blood vessel/smooth muscle",
                             "blood-brain barrier",
                             "blood vessel/smooth muscle",
                             "blood/blood serum",
                             "bone/blood vessel")
    
    KE_annotated$organ_tissue[KE_annotated$organ_tissue %in% serum_related_organs] = "serum"
    
    if(enlarge_ke_to_uncategorized){
      KE_annotated[is.na(KE_annotated$organ_tissue),"organ_tissue"] = tissue_of_interest
    }
    
    unique_tp <- sort(unique(KE_annotated$time)) # time points present in KE enrichment dataframe
    # KE enrichment for each time point in vitro
    
    KE_list_by_tp = list()
    for(tp in unique_tp){
      KE_list_by_tp[[as.character(tp)]] = KE_annotated[KE_annotated$time == tp,]
    }
    
    # --- Molecular weight for current chemical
    chemical_mw <- PFAS_mw[[ss]]
    if (is.null(chemical_mw) || !is.finite(chemical_mw)) {
      stop(sprintf("No molecular weight found for PFAS '%s'.", ss))
    }
    
    # --- Validate Simulation_results structure for this PFAS and tissue
    if (!is.list(Simulation_results[[ss]]) ||
        is.null(Simulation_results[[ss]][[tissue_of_interest]])) {
      stop(sprintf("`Simulation_results[['%s']][['%s']]` is missing.",
                   ss, tissue_of_interest))
    }
    
    
    # avg_pbpk_mod_pred is a dataframe with nrow as the years of PBK simulation and columns as time 
    # of simulation and columns as simulation for each individual PBK model
    # plus one column for each draw of the BMA consensus PBK model
    avg_pbpk_mod_pred =  Simulation_results[[ss]][[tissue_of_interest]]
  
    if (!is.data.frame(avg_pbpk_mod_pred)) {
      stop(sprintf("`Simulation_results[['%s']][['%s']]` must be a data.frame.",
                   ss, tissue_of_interest))
    }
    if (nrow(avg_pbpk_mod_pred) < 1L || ncol(avg_pbpk_mod_pred) < 2L) {
      stop(sprintf("`Simulation_results[['%s']][['%s']]` must have >=1 row and >=2 columns (index + simulations).",
                   ss, tissue_of_interest))
    }
    sim_cols <- avg_pbpk_mod_pred[, -1, drop = FALSE]
    if (!all(vapply(sim_cols, is.numeric, logical(1)))) {
      stop(sprintf("Columns 2+ in `Simulation_results[['%s']][['%s']]` must be numeric.",
                   ss, tissue_of_interest))
    }
    if (any(vapply(sim_cols, function(x) all(is.na(x)), logical(1)))) {
      stop(sprintf("`Simulation_results[['%s']][['%s']]` contains all-NA simulation columns.",
                   ss, tissue_of_interest))
    }
    # `time_point_invivo` as a row index sanity check
    if (time_point_invivo < 1L || time_point_invivo > nrow(avg_pbpk_mod_pred)) {
      stop(sprintf("`time_point_invivo` (%d) is out of bounds for `Simulation_results[['%s']][['%s']]` with %d rows.",
                   time_point_invivo, ss, tissue_of_interest, nrow(avg_pbpk_mod_pred)))
    }
    
    res_triggered_ke = identify_kes_triggered_by_exposure_simulation(PFAS = ss,
                                                                     KE_list_by_tp = KE_list_by_tp,
                                                                     PBPK_predictions=avg_pbpk_mod_pred[,-1],
                                                                     invitro_time =time_point_invitro,
                                                                     tissue_of_interest = tissue_of_interest,
                                                                     chemical_mw = chemical_mw,
                                                                     scaled_time_idx = time_point_invivo,
                                                                     pod_variable = pod_variable)
    
    print(paste("Quantile distributino of nr of triggered KE across PBK samples, for ", PFAS, ", at in vitro time point ", time_point_invitro, ", after ", time_point_invivo, " years PBK simulation", sep = ""))
    print(quantile(unlist(lapply(res_triggered_ke$triggered_KE_names, function(L) length(unique(L))))))
    
    
    # triggered_KE is a list. Length of the list is eaqual to the number of 
    # individual PBK models + number of drawn from BMA models.
    # Each position of the list contains a dataframe of the same format of the 
    # result of KE enrichment, but it only contains the KEs that are triggered
    # when considering the internal concentraition of thae specific PBK result 
    # (either the single PBK model or one of the draw of the BMA model)
    triggered_KE = res_triggered_ke$triggered_KE
    
    
    triggered_ke_list[[ss]] <- triggered_KE
  }
  
  return(list("triggered_ke_list"=triggered_ke_list))
  
}

#' Identify Key Events Triggered by Exposure Simulation
#'
#' @description
#' Given a PFAS identifier, a list of KE tables indexed by in vitro time,
#' a matrix/data frame of PBPK-predicted tissue concentrations over time,
#' and a selected point-of-departure (POD) variable (e.g., BMD), this
#' function determines which Key Events (KEs) are considered triggered at
#' a specific (scaled) time index by comparing concentration thresholds
#' derived from the KE table with the simulated concentration at that time.
#'
#' @details
#' The function:
#' 1. Selects the KE table for the provided in vitro time (`invitro_time`)
#'    from `KE_list_by_tp` (i.e., `KE_list_by_tp[[as.character(invitro_time)]]`).
#' 2. Iterates over PBPK prediction columns .
#' 3. For each PBPK prediction column, retrieves the simulated concentration
#'    at `scaled_time_idx` and filters KEs whose threshold
#'    `as.numeric(kes[, pod_variable]) * chemical_mw` is **less than or equal to**
#'    that concentration.
#' 4. Stores both the full KE rows and just the KE descriptions under names
#'    that encode PFAS, in vitro time, and the PBPK column name.
#'
#' **Assumptions & Requirements**
#' - `KE_list_by_tp` is a named list where names are character representations
#'   of time points (e.g., `"10"`) and each element is a data frame containing
#'   at least the columns given in the example structure, including `Ke_description`
#'   and the POD column specified in `pod_variable` (e.g., `"BMD"`, `"BMDL"`, `"BMDU"`).
#' - `PBPK_predictions` is a data.frame or matrix with rows corresponding to time
#'   points (or simulation indices) and columns corresponding to different
#'   simulation outputs or draws. The first column may be time or any reference
#'   column and is **skipped** by design.
#' - `scaled_time_idx` indexes the **row** of `PBPK_predictions` to use.
#' - `chemical_mw` (molecular weight) is used to scale the POD threshold as
#'   `threshold = POD * chemical_mw`. Ensure that the units of `POD` and the
#'   PBPK concentration are compatible after multiplying by `chemical_mw`.
#' - `tissue_of_interest` is currently **not used** (the line filtering by
#'   `organ_tissue` is commented out); keep for forward compatibility or
#'   re-enable if desired.
#'
#' **Units**
#' - The comparison uses: `(POD * chemical_mw) <= concentration_at_tp_i`.
#'   Ensure that this produces comparable units. If your POD is already
#'   expressed in mass/volume compatible with `PBPK_predictions`, you may not
#'   need to multiply by `chemical_mw`.  
#'
#' @param PFAS Character scalar. The PFAS identifier (e.g., `"PFOS"`), used in
#'   the naming of outputs.
#' @param KE_list_by_tp A named list of data frames, each representing KEs at a
#'   given in vitro time point. The list is keyed by `as.character(invitro_time)`.
#'   Each data frame must include at least the column specified by `pod_variable`
#'   and `Ke_description`. Other columns may include metadata (e.g., `TermID`,
#'   `organ_tissue`, `level`, etc.).
#' @param PBPK_predictions A data frame or matrix of PBPK-derived tissue
#'   concentrations (rows = time indices, columns = simulations/variants).
#'   Column 1 is presumed non-comparable to thresholds and is skipped; columns
#'   2..n are evaluated.
#' @param invitro_time Numeric or character scalar that selects the KE table
#'   from `KE_list_by_tp` via `as.character(invitro_time)`.
#' @param tissue_of_interest Character scalar indicating the tissue of interest
#'   (e.g., `"liver"`). Currently not used; retained for forward compatibility.
#' @param chemical_mw Numeric scalar. Molecular weight of the chemical (e.g., PFAS),
#'   used to scale the POD to a comparable concentration threshold.
#' @param scaled_time_idx Integer scalar. Row index into `PBPK_predictions`
#'   indicating which simulation time point to evaluate.
#' @param pod_variable Character scalar naming the column in the KE table to be
#'   used as the POD (e.g., `"BMD"`, `"BMDL"`, `"BMDU"`). Values must be numeric
#'   or coercible via `as.numeric()`.
#'
#' @return
#' A named list with two elements:
#' \itemize{
#'   \item \code{triggered_KE}: A named list of data frames (one per PBPK
#'   prediction column from 2..n) containing the subset of KEs whose thresholds
#'   are \eqn{\le} the simulated concentration at \code{scaled_time_idx}.
#'   Names follow the pattern \code{<PFAS>_invitro_tp_<invitro_time>_<PBPKcol>}.
#'   \item \code{triggered_KE_names}: A named list of character vectors for
#'   the corresponding \code{Ke_description} values in the same order/naming.
#' }
#'
#' @export
identify_kes_triggered_by_exposure_simulation = function(PFAS,
                                                         KE_list_by_tp,
                                                         PBPK_predictions,
                                                         invitro_time,
                                                         tissue_of_interest,
                                                         chemical_mw,
                                                         scaled_time_idx,
                                                         pod_variable){
  kes = KE_list_by_tp[[as.character(invitro_time)]]

  kes[,pod_variable] = as.numeric(kes[,pod_variable])* chemical_mw
  PODs = kes[,pod_variable]
  
  triggered_KE = vector("list", ncol(PBPK_predictions) )
  triggered_KE_names = vector("list", ncol(PBPK_predictions) )
  
  for(j in 1:ncol(PBPK_predictions)){
    concentration_at_tp_i = PBPK_predictions[scaled_time_idx, j]
    kej = kes[which(PODs <= concentration_at_tp_i),]
    
    triggered_KE[[j]] = kej
    triggered_KE_names[[j]] = kej$Ke_description 
  }
  
  names(triggered_KE) = paste(PFAS,"_invitro_tp_",invitro_time,"_",colnames(PBPK_predictions),sep="" )
  names(triggered_KE_names) = paste(PFAS,"_invitro_tp_",invitro_time,"_",colnames(PBPK_predictions),sep="" )
  
  return(list("triggered_KE"=triggered_KE,"triggered_KE_names"=triggered_KE_names))
}

# optimal_models_stats is the output of bmdx tool
# chemical and experiment are the names of the chemicals and the name of experiments that should be present in optimal_model_stats and Enrichment_data
# Enrichment_data is the result of KE enrichment analysis
# it gives in output a dataframe named pods with genes (column Feature) and their points of departures columns (BMD, BMDL and BMDU)

#' Extract Points of Departure (PODs) by Chemical and Experiment
#'
#' @description
#' Filters BMD modeling results to the genes implicated by a specific
#' KE enrichment result for a given chemical and experiment, returning a
#' PODs table with \code{Feature}, \code{BMD}, \code{BMDL}, and \code{BMDU}.
#'
#' @section Required columns:
#' \describe{
#'   \item{optimal_models_stats}{Must contain columns: \code{Experiment}, \code{time},
#'         \code{Feature}, \code{BMD}, \code{BMDL}, \code{BMDU}.}
#'   \item{Enrichment_data}{Must contain columns: \code{Experiment}, \code{Genes}
#'         (semicolon-separated gene symbols).}
#' }
#'
#' @param optimal_models_stats A data frame from the BMDx tool containing model
#'   statistics and PODs. Must include the columns listed under *Required columns*.
#' @param chemical Character scalar. Name (or pattern) of the chemical used to
#'   subset \code{optimal_models_stats$Experiment} via \code{grep}.
#' @param experiment Character scalar. Exact experiment identifier expected to match
#'   \verb{paste(Experiment, time, sep = "_")} in \code{optimal_models_stats}, and to
#'   match \code{Enrichment_data$Experiment}.
#' @param Enrichment_data A data frame with KE enrichment results, including an
#'   \code{Experiment} column and a semicolon-separated \code{Genes} column.
#'
#' @return
#' A data frame named \code{pods} with columns \code{Feature}, \code{BMD},
#' \code{BMDL}, and \code{BMDU}. Row names are set to \code{Feature}, and rows
#' are ordered to follow the gene order parsed from \code{Enrichment_data$Genes}.
#'
extract_pods_by_chemical_and_experiment = function(optimal_models_stats,chemical, experiment,Enrichment_data ){
  OMS = optimal_models_stats[grep(pattern = chemical,x = optimal_models_stats$Experiment),]

  Enrichment_data = Enrichment_data[Enrichment_data$Experiment %in% experiment,]
  genes = unique(unlist(strsplit(Enrichment_data$Genes,split = ";")))
  OMSexpi = OMS[paste(OMS$Experiment, OMS$time, sep = "_")==experiment,]

  # get all the genes for a specific experiment
  pods = OMSexpi[OMSexpi$Feature %in% genes,c("Feature","BMD","BMDL","BMDU")]
  rownames(pods) = pods$Feature
  pods = pods[genes,]
  return(pods)
}

# this function takes in input 
# NPerm: number of permutations
# pods dataframe with genes (column Feature) and their points of departures columns (BMD, BMDL and BMDU)
# For each gene in the dataframe (each row of the dataframe) the function computes a random BMD in the range BMDL - BMDU by uniform sampling
# It gives  in ouput the matrix permutated_BMD_by_gene which genes on the rows and NPerm + 1 colums with the first column being the original BMD of the genes and the other columns being the randomly sampled BMD

#' Uniform Sampling of BMD within Gene-Specific POD Intervals
#'
#' @description
#' For each gene (row) in a PODs data frame, draw \code{NPerm} random
#' Benchmark Dose (BMD) values uniformly within the interval
#' \code{[BMDL, BMDU]}. Returns a matrix/data frame with genes as rows and
#' \code{NPerm + 1} columns: the first column holds the original \code{BMD}
#' and the remaining columns contain the sampled BMDs.
#'
#' @param NPerm Integer (non-negative). Number of uniform samples to draw per gene.
#' @param pods A data frame with at least the columns \code{Feature}, \code{BMD},
#'   \code{BMDL}, and \code{BMDU}. Each row corresponds to one gene.
#'
#' @return
#' A data frame with \code{nrow(pods)} rows (one per gene) and \code{NPerm + 1}
#' columns:
#' \itemize{
#'   \item \code{originalBMD}: the original \code{pods$BMD}.
#'   \item \code{Perm_1_BMD}, \code{Perm_2_BMD}, \dots, \code{Perm_NPerm_BMD}:
#'   uniform samples within \code{[BMDL, BMDU]} for each gene.
#' }
#' Row names are set to \code{pods$Feature}.
#'
uniform_sampling_BMD_for_KE = function(NPerm, pods, bmd_upper_level="BMDU"){
  genes = pods$Feature

  permutated_BMD_by_gene = matrix(0, nrow = length(genes), ncol = NPerm)
  rownames(permutated_BMD_by_gene) = genes

  for(np in 1:NPerm){
    rgenes = apply(X = pods, MARGIN = 1, FUN = function(riga){
      runif(n = 1, min = as.numeric(riga["BMDL"]), max = as.numeric(riga[bmd_upper_level]))
    })
    permutated_BMD_by_gene[,np]= rgenes
  }

  permutated_BMD_by_gene = cbind("originalBMD" = pods$BMD, permutated_BMD_by_gene)
  colnames(permutated_BMD_by_gene) = c("originalBMD", paste("Perm_",1:NPerm,"_BMD",sep = ""))
  permutated_BMD_by_gene = as.data.frame(permutated_BMD_by_gene)

  return(permutated_BMD_by_gene)
}

# It takes in input the original KE enrichment list (Enrichment_data) in the output format of AOPfingerpint
# It takes in input the matrix permutated_BMD_by_gene which genes on the rows and NPerm + 1 colums with the first column being the original BMD of the genes and the other columns being the randomly sampled BMD
# percentile is the percentile of BMD across the genes in a KE that is used to summarize the BMD from the gene level to KE. Default is 5-th percentile
# it gives in output EnrichmentResult_permutation_list a list of dataframe of length NPerm, each dataframe contains the KE enrichment result with permutated BMD values

#' Create KE Enrichment Lists Using Randomly Sampled Gene-Level BMD Values
#'
#' @description
#' For each permutation (including the original), summarize gene-level BMD values
#' into KE-level BMDs by taking a chosen percentile across genes in each KE.
#' Returns a list of KE enrichment data frames—one per permutation—where the
#' \code{BMD} column has been replaced by the percentile-based summary computed
#' from \code{permutated_BMD_by_gene}.
#'
#'
#' @param Enrichment_data A KE enrichment data frame (AOPfingerprint output) that
#'   includes \code{Experiment} and \code{Genes} columns. May contain additional
#'   metadata columns that will be carried through to the output.
#' @param permutated_BMD_by_gene A matrix or data frame with row names equal to
#'   gene identifiers and \code{NPerm + 1} columns where the first column holds
#'   the original BMD values and the remaining columns hold permuted BMD values.
#' @param NPerm Integer (non-negative). Number of permutation columns present
#'   beyond the original BMD column.
#' @param experiment Character scalar. Experiment identifier used to subset
#'   \code{Enrichment_data}.
#' @param percentile Numeric in \code{[0, 1]} (default \code{0.05}). The
#'   percentile of gene-level BMD values within each KE that summarizes to a
#'   KE-level BMD (e.g., 5th percentile).
#'
#' @return
#' A list of length \code{NPerm + 1}. Each element is a data frame mirroring
#' \code{Enrichment_data[Enrichment_data$Experiment \%in\% experiment, ]}, but with
#' the \code{BMD} column replaced by the percentile summary of the corresponding
#' column in \code{permutated_BMD_by_gene}:
#' \itemize{
#'   \item List element \code{[[1]]}: KE BMDs computed from \code{originalBMD}.
#'   \item List elements \code{[[2]]}..\code{[[NPerm+1]]}: KE BMDs computed from
#'         permuted columns.
#' }
#'
create_enrichment_list_with_randomly_sampled_BMD_values = function(Enrichment_data, 
                                                                             permutated_BMD_by_gene, 
                                                                             NPerm, 
                                                                             experiment,
                                                                             percentile = 0.05){
  Enrichment_data = Enrichment_data[Enrichment_data$Experiment %in% experiment,]
  
  # Pre-split Genes ONCE  
  gene_list <- strsplit(Enrichment_data$Genes, ";", fixed = TRUE)
  
  # Pre-map gene symbols -> row indices ONCE (preserves duplicates + order)
  rn <- rownames(permutated_BMD_by_gene)
  idx_lookup <- setNames(seq_along(rn), rn)
  
  gene_idx_list <- lapply(gene_list, function(gi) unname(idx_lookup[gi]))
  
  # Preallocate output list
  EnrichmentResult_permutation_list <- vector("list", NPerm + 1L)
  
  for (np in 1:(NPerm + 1)) {
    EDi <- Enrichment_data
    
    bmd_col <- permutated_BMD_by_gene[, np]
    
    # Compute BMD for each KE 
    for (ke in 1:nrow(Enrichment_data)) {
      gi <- gene_idx_list[[ke]]
      th_05_BMD <- quantile(bmd_col[gi], percentile)  
      EDi[ke, "BMD"] <- th_05_BMD
    }
    EnrichmentResult_permutation_list[[np]] <- EDi
  }
  return(EnrichmentResult_permutation_list)
}


# This function takes in input:
# pfas: a string that can be either PFOS or PFOA
# exposure_time: an integer, representing the number of years of PBK simulation
# tissue_of_interest: a string representing the PBK compartmnet of interest (e.g. liver, serum)
# pod_variable: a stringthat represent the point of departure value (e.g. BMD)
# time_point_invitro: an integer representing the time point of the invitro data
# Simulation_results: nested list containing the PBK simulation results
# list_enrichment: a list of KE enrichment dataframes, one for each PFAS
# enlarge_ke_to_uncategorized: boolean. If TRUE, KE with no specific organ of relevance assigned are also included in the analysis
# it gives in output a 3D tensor T of dimensions: k- key events, d-draw from BMA 
# model, t time of PBK simulation. Each position is T[k, d, t] = BMD value for friggered KE k, draw d, at time t


#' Find KE Activation Over PBK Simulation Time
#'
#' For a given PFAS (e.g., \code{"PFOS"} or \code{"PFOA"}), tissue, and point-of-departure
#' variable (e.g., BMD), this function integrates PBPK simulation outputs with omics-based
#' enrichment results (via \code{omics_pbpk_integration()}) at each PBK time point and
#' constructs a 3D tensor with BMD values of triggered Key Events (KEs) across BMA draws
#' and PBK time.
#'
#' Specifically, at each time point \code{t = 1, ..., exposure_time}, it extracts—per BMA draw—
#' the set of triggered KEs and their BMD values, arranges them into a KE × draw matrix
#' (filling non-triggered KEs with 0), then pads across time to the union of KEs observed
#' at any time point. Finally, it stacks these matrices into a 3D array
#' \code{[KE × draw × time]}.
#'
#' @param pfas Character scalar. The PFAS identifier, typically \code{"PFOS"} or \code{"PFOA"}.
#' @param exposure_time Integer scalar. Number of PBK simulation time points (e.g., years).
#' @param tissue_of_interest Character scalar. PBK compartment of interest (e.g., \code{"liver"}, \code{"serum"}).
#' @param pod_variable Character scalar. Point-of-departure variable to extract (e.g., \code{"BMD"}).
#' @param time_point_invitro Integer scalar. Time point of the in vitro data used in the integration.
#' @param Simulation_results Nested list with PBK simulation outputs required by
#'   \code{omics_pbpk_integration()}.
#' @param list_enrichment A list of KE enrichment data frames, indexed by PFAS name,
#'   to be supplied to \code{omics_pbpk_integration()} as \code{Enrichment_data}.
#' @param enlarge_ke_to_uncategorized Logical; if \code{TRUE}, include KEs without a specific
#'   organ/tissue label in the analysis. Default is \code{FALSE}.
#'
#' @return A 3D numeric array (tensor) with dimensions:
#'   \itemize{
#'     \item \strong{Dimension 1 (rows)}: KE names (union across all time points).
#'     \item \strong{Dimension 2 (columns)}: BMA draw identifiers (e.g., \code{"draw_1"}, \code{"draw_2"}, ...).
#'     \item \strong{Dimension 3 (slices)}: PBK time points \code{1..exposure_time}.
#'   }
#'   Each entry \code{T[ke, draw, t]} is the BMD value for KE \code{ke} triggered by BMA \code{draw}
#'   at PBK time \code{t}; non-triggered KEs are represented by \code{0}.
#'
find_ke_activation_over_time = function(pfas,
                                      exposure_time,
                                      tissue_of_interest,
                                      pod_variable,
                                      time_point_invitro,
                                      Simulation_results,
                                      list_enrichment,
                                      enlarge_ke_to_uncategorized = F){
ke_by_tp = vector("list", length(exposure_time))

for(tpi in 1:exposure_time){
    
    # PFAS: a string that can be either PFOS or PFOA
    # tissue_of_interest: a string representing the PBK compartmnet of interest (e.g. liver, serum)
    # pod_variable: a stringthat represent the point of departure value (e.g. BMD)
    # time_point_invitro: an integer representing the time point of the invitro data
    # tpi: integer, one of the exposure time of PBK simulation
    # Simulation_results: nested list containing the PBK simulation results
    # Enrichment_data: a list of KE enrichment dataframes, one for each PFAS
    # enlarge_ke_to_uncategorized: if true, it includes also KE that do not have a tissue label
    
    res = omics_pbpk_integration(PFAS=pfas,
                                 tissue_of_interest=tissue_of_interest,
                                 pod_variable=pod_variable,
                                 time_point_invitro=time_point_invitro,
                                 time_point_invivo=tpi,
                                 Simulation_results=Simulation_results,
                                 Enrichment_data=list_enrichment, #the np-th permutated BMD
                                 enlarge_ke_to_uncategorized = enlarge_ke_to_uncategorized)
    
    pfas_triggered_ke_list = res$triggered_ke_list[[pfas]]
    n_draw = length(grep(pattern = "draw_",x = names(pfas_triggered_ke_list))) # Nr of drawn in the BMA algorithm
    
    # bmdi_list is a list where for each draw from the BMD model we store the 
    # BMD of the KEs that are triggered by the comparison with the internal 
    # concentration predicted by draw_i sample
    bmdi_list <- vector("list", length = n_draw)
    for(draw_i in 1:n_draw){
      ddi = unique(pfas_triggered_ke_list[[paste(pfas,"_invitro_tp_",time_point_invitro,"_C",tissue_of_interest,"_BMA_draw_",draw_i,sep = "")]][,c("Ke_description","BMD")])
      ketpi = ddi$Ke_description
      bmdi = as.numeric(ddi$BMD)
      names(bmdi) = ketpi
      bmdi_list[[draw_i]] = bmdi
    }
    
    # all_names is the list of al KEs triggered by at list one draw of the BMA model
    all_names <- unique(unlist(lapply(bmdi_list, names)))
    
    # bmdi_mat is a matrix with KEs on the rows and the nr of BMA draw on the colums
    # for each KE and draw, the matrix contain the BMD value of the KEs triggered
    # by the draw sample of the BMA model
    bmdi_mat <- matrix(0, nrow = length(all_names), ncol = length(bmdi_list),
                       dimnames = list(all_names, paste0("draw_", seq_along(bmdi_list))))
    
    for (i in seq_along(bmdi_list)) {
      vec <- bmdi_list[[i]]
      bmdi_mat[names(vec), i] <- vec
    }
    
    ke_by_tp[[tpi]] = bmdi_mat
  }
  
  # ke_by_tp list is long as the years of PBK simulation.
  # for each year, it contains the bmdi mat.
  # bmdi_mat is a matrix with KEs on the rows and the nr of BMA draw on the colums
  # for each KE and draw, the matrix contain the BMD value of the KEs triggered
  # by the draw sample of the BMA model
  
  # all_rows <- unique(unlist(lapply(ke_by_tp, rownames)))
  # 
  # # if no KE pass the filters, all_rows will be NULL. 
  # if(is.null(all_rows)){
  #   all_rows = unique(list_enrichment[[pfas]]$TermID)
  # }
  all_rows = unique(list_enrichment[[pfas]]$Ke_description)
  all_cols <- colnames(ke_by_tp[[1]])
  stopifnot(all(sapply(ke_by_tp, function(m) all(colnames(m) == all_cols))))
  
  ke_by_tp_padded <- lapply(ke_by_tp, function(m) {
    tmp <- matrix(0, nrow = length(all_rows), ncol = length(all_cols),
                  dimnames = list(all_rows, all_cols))
    tmp[rownames(m), ] <- m
    tmp
  })
  
  xs = unlist(lapply(ke_by_tp_padded, nrow))
  print("Nr of KEs in ke_by_tp_padded--->")
  print(xs)
  print(all(xs == xs[1]))
  print("Nr of KE actually triggered")
  print(length(unique(unlist(lapply(ke_by_tp, rownames)))))
  
  ke_tensor <- simplify2array(ke_by_tp_padded)
  
  # ke_tensor is a 3D tensor, with first dimension, nr of KEs, second dimension
  # nr of BMA draw samples and 3rd dimension, nr of PBK year simulation
  # Each position of the tensor is the BMD value of the KE
  return(ke_tensor)
}


# given the list of KEs with permuted BMD values assigned, this function compute 
# the probability of KE activation for each one of its permuted BMD given the
# internal concentration of the compartment computed by the PBK model
# it takes in input EnrichmentResult_permutation_list a list of dataframe of length NPerm, 
# each dataframe contains the KE enrichment result with permutated BMD values

#' Compute KE Activation Probability over Permuted BMD Sets
#'
#' @description
#' For each permutation of KE-level BMD values, compute the probability of KE
#' activation over time given internal concentrations from a PBK simulation.
#' Iterates across \code{NPerm} permutations contained in
#' \code{EnrichmentResult_permutation_list} and calls
#' \code{find_ke_activation_over_time()} to produce one result per permutation.
#'
#' @param EnrichmentResult_permutation_list A list of KE enrichment data frames,
#'   each corresponding to one permutation and containing KE-level \code{BMD}
#'   values to be evaluated for activation.
#' @param pfas Character scalar. Identifier (e.g., compound name) used as the
#'   list key for \code{list_enrichment} passed into
#'   \code{find_ke_activation_over_time()}.
#' @param NPerm Integer (positive). Number of permutations to iterate over. The
#'   function processes list elements \code{[[1]]} through \code{[[NPerm]]}.
#' @param exposure_time Numeric or character. Exposure duration or schedule
#'   required by \code{find_ke_activation_over_time()}.
#' @param tissue_of_interest Character scalar. Tissue/compartment in which KE
#'   activation is evaluated (must be consistent with \code{Simulation_results}).
#' @param pod_variable Character scalar. Name of the POD variable used inside
#'   \code{find_ke_activation_over_time()} (e.g., \code{"BMD"}).
#' @param time_point_invitro Numeric or character. In vitro reference time point
#'   (or mapping) used by \code{find_ke_activation_over_time()}.
#' @param Simulation_results PBK simulation outputs with internal concentrations
#'   over time, formatted as required by \code{find_ke_activation_over_time()}.
#' @param enlarge_ke_to_uncategorized Logical (default \code{TRUE}). Passed through
#'   to \code{find_ke_activation_over_time()} to optionally include
#'   "uncategorized" KEs.  
#'
#' @return
#' A list \code{M_list} of length \code{NPerm}, where each element is the return
#' value of \code{find_ke_activation_over_time()} for the corresponding
#' permutation. The exact structure of each element depends on
#' \code{find_ke_activation_over_time()}, typically containing KE activation
#' probabilities over time for the requested tissue and exposure scenario.
#'
KE_activation_over_permutated_BMD = function(EnrichmentResult_permutation_list, 
                                             pfas, 
                                             NPerm,
                                             exposure_time,
                                             tissue_of_interest,
                                             pod_variable,
                                             time_point_invitro,
                                             Simulation_results,
                                             enlarge_ke_to_uncategorized = T){
  
  M_list = vector("list", NPerm)
  pb = txtProgressBar(min = 1, max = NPerm, style = 3)
  
  for(np in 1:NPerm){
    print(paste("Permutation nr:", np))
    list_enrichment = list()
    list_enrichment[[pfas]] = EnrichmentResult_permutation_list[[np]]
    M = find_ke_activation_over_time(pfas,
                                     exposure_time,
                                     tissue_of_interest,
                                     pod_variable,
                                     time_point_invitro,
                                     Simulation_results,
                                     list_enrichment,
                                     enlarge_ke_to_uncategorized = T)
    
    print(class(M))
    
    M_list[[np]] = M
    
    print(class( M_list[[np]]))
    
    setTxtProgressBar(pb, np)
  }
  close(pb)
  
  return(M_list)
}

# given the M_list object from function KE_activation_over_permutated_BMD, this 
# function computes the probability of activation of each KE across the BMA
# draw sampels, across the NPerm BMD permutation rounds


#' Compute KE Activation Percentages Across BMD Permutations and BMA Samples
#'
#' Given a list of 3D arrays produced by `KE_activation_over_permutated_BMD()`,
#' this function computes, for each Key Event (KE), the percentage of activation
#' across the list elements (e.g., permutations or samples) while preserving the
#' original array dimensions (KE × BMA samples × BMD permutation rounds).
#'
#' The function first aligns KE rows across all arrays by padding any missing KEs
#' with zeros (non-activated), using the KE universe provided in
#' `Biological_system_annotations$key_event_name`. It then treats any value `> 0`
#' as "activated" and calculates, for each array cell, the percentage of list
#' elements where activation occurs.
#'
#' @param M_list A list of 3D numeric arrays. Each array must have:
#'   \itemize{
#'     \item \strong{Dimension 1 (rows)} = KEs, with rownames set to KE names.
#'     \item \strong{Dimension 2} = typically BMA draw samples (must be compatible across list elements).
#'     \item \strong{Dimension 3} = typically BMD permutation rounds (must be compatible across list elements).
#'   }
#'   All arrays in the list must share the same dimensions for the 2nd and 3rd axes
#'   and (ideally) identical \code{dimnames} for these axes. The 1st axis (KEs)
#'   is harmonized internally by padding missing KEs with zeros according to the KE
#'   universe.
#' @param NPerm Number of permutation
#'
#' @return A 3D numeric array of the same shape as the padded inputs:
#'   \code{[length(Biological_system_annotations$key_event_name) × dim2 × dim3]},
#'   where each element is a percentage in \code{[0, 100]} indicating the proportion
#'   (over the list elements of \code{M_list}) where activation (\code{> 0}) occurs
#'   at that specific KE × sample × permutation index. \code{dimnames} are preserved
#'   where available and re-ordered to match the KE universe.
#'
#' @details
#' \itemize{
#'   \item \strong{KE universe:} The set of all KEs is taken from
#'   \code{Biological_system_annotations$key_event_name}, which must be a character
#'   vector available in the current environment. Any KE missing in a given array
#'   is padded with a zero-plane along the first dimension.
#'
#'   \item \strong{Activation rule:} Values \code{> 0} are considered "activated".
#'   Values \code{\link[base]{NA}} are ignored in the counting step (i.e., they do not
#'   increment the activation count). For robust behavior, it is recommended to
#'   ensure non-activated entries are explicitly set to \code{0}.
#'
#'   \item \strong{Across-list aggregation:} For each cell \code{[i, j, k]},
#'   the function counts how many list elements have \code{> 0} and divides by
#'   \code{length(M_list)}, returning a percentage.
#'
#'   \item \strong{Dimensional consistency:} The function assumes that the 2nd and 3rd
#'   dimensions (and their \code{dimnames}) are consistent across all arrays. If they
#'   differ, results may be misaligned or \code{dimnames} may not reflect a true union.
#' }
#'
#' @section Input Structure:
#' Each element \code{M} in \code{M_list} is expected to be a 3D numeric array:
#' \preformatted{
#'   dim(M)      = c(n_ke_M, n_samples, n_permutations)
#'   dimnames(M) = list(KE_names, sample_names, permutation_names)
#' }
#' The KE names (\code{rownames}) are aligned to
#' \code{Biological_system_annotations$key_event_name}.
#'
#' @importFrom abind abind
compute_percentages_of_KE_activation_across_BMD_permutation_and_BMA_samples = function(M_list, NPerm){
  
  # The function first check that M_list includes all KEs for each permutation.
  # If not, pads the M_list matrices to have all KEs and add zeros for the KE 
  # that were not present
  
  # Get the union of all KE (rownames) across all matrices
  all_kes <- unique(unlist(lapply(M_list, rownames)))
  # all_kes = Biological_system_annotations$key_event_name
  
  print("str(M_list)")
  print(str(M_list))
  # Pad each matrix to include all KEs, filling missing ones with zeros
  M_list_padded <- lapply(M_list, function(M) {
    
    # Determine missing KEs in this matrix
    missing_kes <- setdiff(all_kes, rownames(M))
    # Create zero array for missing KEs
    if (length(missing_kes) > 0) {
      zero_block <- array(0, dim = c(length(missing_kes), dim(M)[2], dim(M)[3]))
      dimnames(zero_block) <- list(missing_kes,
                                   dimnames(M)[[2]],
                                   dimnames(M)[[3]])
      
      # Combine existing and missing
      M_full <- abind::abind(M, zero_block, along = 1)
    } else {
      M_full <- M
    }
    
    # Reorder rows to match all_kes
    M_full <- M_full[all_kes, , , drop = FALSE]
    
    return(M_full)
  })
  
  # After, the probability of activation of each KE is computed
  
  M_count = M_list_padded[[1]]
  M_count[M_count>0] = 1
  
  if(length(M_list_padded)>1){
    for(index in 2:length(M_list_padded)){
      mi = M_list_padded[[index]]
      M_count[mi>0] = M_count[mi>0] + 1
    }
  }
  
  # M_count = ( M_count/length(M_list_padded) ) * 100
  M_count = ( M_count/NPerm ) * 100
  
}

#' Compute PBK–BMD KE Activation Probability Matrices
#'
#' For each tissue and in vitro experiment, this function estimates the probability
#' (in percent) that a Key Event (KE) is activated across \code{NPerm} resamplings
#' of gene-level BMD within their \code{[BMDL, BMDU]} ranges. The workflow is:
#' \enumerate{
#'   \item Extract experiment-specific POD/BMD statistics.
#'   \item Perform \code{NPerm} uniform resamplings of gene-level BMD between
#'         \code{BMDL} and \code{BMDU}.
#'   \item Recompute KE-level enrichment for each permutation, summarizing genes'
#'         BMDs to a KE-level statistic via the provided percentile.
#'   \item Run PBPK/omics integration across PBK time points and BMA draws to obtain,
#'         per permutation, KE activation over time.
#'   \item Collapse across permutations to probabilities (\%), resulting in a
#'         \code{KE × draw × time} array per tissue/experiment.
#' }
#'
#' @param tissues_of_interest Character vector of PBK compartments of interest
#'   (e.g., \code{c("liver","serum")}).
#' @param unique_exp Character vector of experiment identifiers present in the in vitro
#'   dataset, typically of the form \code{"CHEMICAL_TIME"} (e.g., \code{"PFOA_14"}).
#' @param optimal_models_stats Result object from BMD analysis on all experiments in the
#'   in vitro data (used by \code{extract_pods_by_chemical_and_experiment()}).
#' @param Enrichment_data Result of KE enrichment analysis based on BMD output (used as
#'   the template for permutation-based enrichment recomputation).
#' @param NPerm Integer; number of uniform resamplings in \code{[BMDL, BMDU]} for genes'
#'   BMD values used to associate genes with KEs.
#' @param BMD_summarization_percentile Numeric; percentile used to summarize gene-level
#'   resampled BMDs to a KE-level value (e.g., \code{5} for the 5th percentile).
#' @param exposure_time Integer; number of PBK simulation time points (e.g., years).
#' @param pod_variable Character scalar specifying the point-of-departure variable used
#'   for activation computation; for the current setup this should be \code{"BMD"}.
#' @param Simulation_results PBK simulation results object required by
#'   \code{KE_activation_over_permutated_BMD()} and downstream integration.
#'
#' @return A nested named list \code{probability_matrices} such that:
#' \preformatted{
#' probability_matrices[[tissue]][[experiment]] -> 3D numeric array
#'   with dimensions [n_ke × n_draw × exposure_time]
#' }
#' where each entry is a percentage in \code{[0, 100]} representing the probability
#' that KE \code{ke} is activated for BMA \code{draw} at PBK time \code{t}, across the
#' \code{NPerm} permutations of gene-level BMD. The inner array inherits \code{dimnames}
#' for KEs and draws where available; time slices correspond to \code{1..exposure_time}.
#'
#' @details
#' \itemize{
#'   \item \strong{Dependencies:} This function orchestrates several helpers:
#'   \code{extract_pods_by_chemical_and_experiment()},
#'   \code{uniform_sampling_BMD_for_KE()},
#'   \code{create_enrichment_list_with_randomly_sampled_BMD_values()},
#'   \code{KE_activation_over_permutated_BMD()}, and
#'   \code{compute_percentages_of_KE_activation_across_BMD_permutation_and_BMA_samples()}.
#'   \item \strong{Activation definition:} Downstream, activation is typically defined as
#'   values \code{> 0} in the KE × draw × time arrays, which are aggregated to percentages
#'   over permutations.
#'   \item \strong{Progress messages:} The function prints the current \code{tissue} and
#'   \code{experiment} being processed for basic progress reporting.
#' }
#'
#' @export
compute_PBK_BMD_probability_matrix_activation = function(tissues_of_interest,
                                                         unique_exp,
                                                         optimal_models_stats,
                                                         Enrichment_data,
                                                         NPerm,
                                                         BMD_summarization_percentile,
                                                         exposure_time,
                                                         pod_variable = "BMD",
                                                         Simulation_results,
                                                         bmd_upper_level = "BMDU"){
  probability_matrices = list()
  
  # Preallocate outer list (one element per tissue), with names
  probability_matrices <- setNames(vector("list", length(tissues_of_interest)), tissues_of_interest)
  
  # For each tissue, preallocate inner list (one element per experiment), with names
  for (t in tissues_of_interest) {
    probability_matrices[[t]] <- setNames(vector("list", length(unique_exp)), unique_exp)
  }
  
  # Precompute permutated enrichment for each experiment ID - these are tissue independent
  EnrichmentResult_permutation_list_by_experiment <- setNames(vector("list", length(unique_exp)), unique_exp)
  
  print("Permute Samples... ")
  for(experiment_i in unique_exp){
    chemical = strsplit(x = experiment_i, split = "_")[[1]][1]
    
    pods = extract_pods_by_chemical_and_experiment(optimal_models_stats = optimal_models_stats,
                                                   chemical = chemical, 
                                                   experiment = experiment_i,
                                                   Enrichment_data =Enrichment_data)
    
    PermutedSamples = uniform_sampling_BMD_for_KE(NPerm, pods, bmd_upper_level=bmd_upper_level)

    
    EnrichmentResult_permutation_list = create_enrichment_list_with_randomly_sampled_BMD_values(Enrichment_data, 
                                                                                                          PermutedSamples, 
                                                                                                          NPerm, 
                                                                                                          experiment_i,
                                                                                                          percentile = BMD_summarization_percentile)
    
    EnrichmentResult_permutation_list_by_experiment[[experiment_i]] = EnrichmentResult_permutation_list
  }
  
  print("Done!")
  
  for(tissue  in tissues_of_interest ){
    for(experiment_i in unique_exp){
      print(tissue)
      print(experiment_i)
      chemical = strsplit(x = experiment_i, split = "_")[[1]][1]
      time_point_invitro = as.numeric(unlist(strsplit(experiment_i, "_"))[2])

      M_list = KE_activation_over_permutated_BMD(EnrichmentResult_permutation_list = EnrichmentResult_permutation_list_by_experiment[[experiment_i]], 
                                                 pfas= chemical, 
                                                 NPerm,
                                                 exposure_time,
                                                 tissue_of_interest = tissue,
                                                 pod_variable,
                                                 time_point_invitro,
                                                 Simulation_results,
                                                 enlarge_ke_to_uncategorized = T)
      
      to_rem = which(unlist(lapply(M_list, is.list)))
      
      print("Nr of permutation with errors")
      print(to_rem)
      
      if(length(to_rem)>0){
        print(paste("Iterations removed: ", length(to_rem)))
        M_list = M_list[-to_rem]
      }
      
      if(length(M_list)==0){
        next()
      }
      
      M_count = compute_percentages_of_KE_activation_across_BMD_permutation_and_BMA_samples(M_list, NPerm)
      
      # M_counts has percentage values of how many times across random sampling the BMD of a KE is below the internal concentration estimated by each sample of the BMA
      probability_matrices[[tissue]][[experiment_i]] = M_count 
    }
  }
  
  return(probability_matrices)
}

#' Convert Probability Tensors into Long Data Frames
#'
#' This function takes a list of 3D probability tensors (one per experimental
#' condition), averages permutations, reshapes them into long-format data
#' frames, attaches PFAS identifiers and time labels, and merges biological
#' tissue annotations.
#'
#' Each tensor in `probability_matrices` is expected to have dimensions
#' *(n_key_events × n_samples × n_timepoints)*. The function computes, for each
#' timepoint, the mean probability across permutations/samples, reshapes the
#' resulting matrix with `reshape2::melt()`, and aggregates all experiments into
#' one tidy data frame.
#'
#' @param probability_matrices A named list of 3D numeric arrays (tensors),
#'   where each name encodes an experiment in the format `"PFAS_time"`.
#'   Each element must be an array with dimensions
#'   \code{[key events × samples × time points]}
#'
#' @param time_points Numeric vector. The ordered time points associated with
#'   the third dimension of the tensors (default: \code{c(1, 4, 10, 14)}).
#'
#' @param time_unit Character string giving the unit label applied to time
#'   points (default: \code{"Day"}). This value becomes the name of the
#'   corresponding time column in the output data frame.
#'
#' @return A data frame in long format containing:
#'   \itemize{
#'     \item \code{KeyEvent} — key event identifier
#'     \item \code{Year} — time index from melt()
#'     \item \code{Prob} — mean probability across permutations
#'     \item \code{PFAS} — PFAS identifier
#'     \item \code{<time_unit>} — time annotation (e.g. `"Day 4"`)
#'     \item \code{organ_tissue} — merged tissue information from
#'       \code{Biological_system_annotations}
#'   }
#'
#' @details
#' The function expects that experiment names in `probability_matrices`
#' follow the format \code{"PFAS_time"}, where `time` corresponds to the
#' index in `time_points`. The function merges key event annotations from
#' the global object \code{Biological_system_annotations}.
#'
#' @import reshape2
#' @import AOPfingerprintR
#'
#' @export
convert_tensor_to_dataframes = function(probability_matrices,
                                            time_points = c(1,4,10,14), 
                                            time_unit = "Day"){#}, tissue = "liver"){

  df = NULL
  unique_exp = names(probability_matrices)
  for(exi in unique_exp){
    m1 = probability_matrices[[exi]]
    # m1 = m1/100

    n_time = dim(m1)[3]
    mean_over_samples = c()
    for(tp in 1:n_time){
      mean_over_samples = cbind(mean_over_samples,rowMeans(m1[,,tp]))
    }
    m1 = mean_over_samples
    #

    pps = strsplit(exi, split = "_")[[1]][1]
    day = strsplit(exi, split = "_")[[1]][2]

    # df1 <- cbind(melt(m1), "pfas"=paste(pps, " d. ", day, sep = ""))
    # df1 <- cbind(melt(m1), "pfas"=pps, print0(time_unit)=paste(time_unit," ", day, sep = ""))
    
    df1 <- cbind(melt(m1), "pfas" = pps)
    df1[ time_unit ] <- paste(time_unit, day, sep = " ")

    df = rbind(df, df1)
  }

  colnames(df) = c("KeyEvent","Year","Prob","PFAS", time_unit)
  df[,time_unit] = factor(df[,time_unit], levels = paste(time_unit," ", sort(time_points,decreasing = F),sep=""))#c("Day 1", "Day 4", "Day 10", "Day 14"))
  # df$Day = factor(df$Day, levels = c("Day 1", "Day 4", "Day 10", "Day 14"))

  Mapping = unique(AOPfingerprintR::aop_ke_table_hure[,c("Ke","Ke_description")])
  rownames(Mapping) = Mapping$Ke_description
  df$Ke = Mapping[as.character(df$KeyEvent),"Ke"]
  
  df = merge(x = df,
             y = unique(AOPfingerprintR::Biological_system_annotations[,c("ke","organ_tissue")]),
             by.x = "Ke",
             by.y = "ke",
             all.x = TRUE)

  # df$organ_tissue[is.na(df$organ_tissue)] = paste("Systemic", tissue,sep = "_")
  return(df)
}

