
#' Batch wrapper for stochastic reverse dosimetry across PFAS, time, and compartment
#'
#' @description
#' Runs a stochastic reverse dosimetry workflow across combinations of PFAS
#' chemicals, time points, and biological compartments. The function filters
#' the input omics enrichment table by PFAS and time, rescales the POD
#' (point-of-departure) by the molecular weight of the PFAS, and calls
#' \code{reverse_dosimetry()} to estimate human equivalent doses (HEDs).
#' It stores both per-sample HED results and summarized (median/variance)
#' results for each combination.
#'
#' @param Enrichment_data_BMD \code{data.frame}. Omics enrichment table containing,
#' at minimum, the columns:
#' \itemize{
#'   \item \code{PFAS} (character): PFAS identifier (e.g., "PFOA", "PFOS").
#'   \item \code{time} (numeric/integer): Time point to subset.
#'   \item \code{Ke_description} (character): Key event description (will be renamed to \code{key_event}).
#'   \item a column named as in \code{pod} (numeric): POD/BMD values to be rescaled by PFAS molecular weight.
#'   \item \code{TermID} (optional): Used only to deduplicate rows; dropped before modeling.
#' }
#'
#' @param duration \code{numeric}. Exposure duration (same unit expected by
#' \code{reverse_dosimetry()}). Default: \code{40}.
#'
#' @param BW \code{numeric}. Body weight used by the reverse dosimetry model (if applicable
#' to your \code{reverse_dosimetry()} implementation). Default: \code{70}.
#'
#' @param PFAS_mw \code{named numeric}. Named vector mapping PFAS identifiers to their
#' molecular weights (e.g., \code{c("PFOA" = 414.07, "PFOS" = 500.13)}). Used to rescale
#' POD values before modeling.
#'
#' @param pod \code{character}. Name of the column in \code{Enrichment_data_BMD} that
#' contains the POD/BMD values to use. Default: \code{"BMD"}.
#'
#' @param PFAS \code{character}. Vector of PFAS identifiers to iterate over. Must match
#' the names used in \code{PFAS_mw} and \code{Enrichment_data_BMD$PFAS}. Default:
#' \code{c("PFOA","PFOS")}.
#'
#' @param time \code{numeric/integer}. Vector of time points to iterate over. Must match
#' values present in \code{Enrichment_data_BMD$time}. Default: \code{c(1, 4, 10, 14)}.
#'
#' @param compartment \code{character}. Vector of compartments to iterate over (e.g.,
#' \code{c("liver","serum")}). These are added to the modeling data but are not used to
#' subset \code{Enrichment_data_BMD}. Default: \code{c("liver","serum")}.
#' 
#' @param  isStochastic = TRUE
#' @param  isParallel = T 
#' @param  Nsamples = 10 Number of samples for the stochastic reverse dosimetry
#'
#' @details
#' For each combination of \code{pfas} in \code{PFAS}, \code{time_point} in \code{time},
#' and \code{compartment}, the function:
#' \enumerate{
#'   \item Subsets \code{Enrichment_data_BMD} by \code{PFAS == pfas} and \code{time == time_point}.
#'   \item Keeps only \code{TermID}, \code{Ke_description}, \code{PFAS}, and the \code{pod} column,
#'         deduplicates rows, renames to \code{key_event} and \code{POD}, coerces \code{POD} to numeric,
#'         and adds \code{time} and \code{compartment} columns.
#'   \item Multiplies \code{POD} by the specified PFAS molecular weight.
#'   \item Calls \code{reverse_dosimetry(..., isStochastic = TRUE, isParallel = TRUE, Nsamples = 10)}.
#'   \item Extracts the \code{bma} results from each stochastic sample:
#'         \itemize{
#'           \item \strong{Per-sample results} (\code{X}): rbind of each sample's \code{bma} table.
#'           \item \strong{Summarized results} (\code{Y}): per combination, medians of
#'                 \code{exposure}, \code{rel_error}, and \code{POD}, plus \code{sd} and \code{var} of
#'                 \code{exposure} across samples.
#'         }
#'   \item Stores both tables in named lists keyed by \code{"PFAS_time_compartment"}.
#' }
#'
#' @return
#' \strong{Note:} As written, the function does not explicitly return a value.
#' For practical use, consider adding:
#' \preformatted{
#'   return(list(
#'     reverse_dosimetry_list = reverse_dosimetry_list,
#'     reverse_dosimetry_list_summarized = reverse_dosimetry_list_summarized
#'   ))
#' }
#'
#' If this return statement is added, the function returns a named list with:
#' \itemize{
#'   \item \code{reverse_dosimetry_list}: a list of per-sample HED result tables (\code{X}).
#'   \item \code{reverse_dosimetry_list_summarized}: a list of summarized HED tables (\code{Y})
#'         containing medians and dispersion statistics per combination.
#' }
#'
#' @section Expected structure of \code{reverse_dosimetry()} output:
#' This wrapper expects \code{reverse_dosimetry()} to return a list accessible with
#' \code{rd[[pfas]]}, where each element is one stochastic sample, and each sample contains
#' a \code{bma} component that is either a data.frame or a list convertible via
#' \code{do.call(cbind, ...)}. The \code{bma} table is expected to include columns
#' \code{exposure}, \code{rel_error}, and \code{POD}.
#'
#' @note
#' \itemize{
#'   \item The call \code{reverse_dosimetry(chemical = pfas, BW = PFAS_mw[pfas], ...)} passes the
#'         PFAS molecular weight to the \code{BW} argument. If \code{BW} is intended to be body weight,
#'         you may want to pass \code{BW = BW} instead.
#'   \item The \code{compartment} values are appended to the modeling data but are not used to subset
#'         \code{Enrichment_data_BMD}. This is by design in the current implementation.
#'   \item The function prints each loop index (\code{pfas}, \code{time_point}, \code{compartment})
#'         for progress tracking.
#' }
#'
#' @examples
#' \dontrun{
#' # Example data (toy)
#' df <- data.frame(
#'   TermID = paste0("T", 1:6),
#'   Ke_description = c("KE1","KE2","KE3","KE1","KE2","KE3"),
#'   PFAS = c("PFOA","PFOA","PFOA","PFOS","PFOS","PFOS"),
#'   time = c(1,1,4,1,4,4),
#'   BMD = c(0.5, 0.8, 1.2, 0.3, 0.6, 0.9)
#' )
#'
#' PFAS_mw <- c("PFOA" = 414.07, "PFOS" = 500.13)
#'
#' # Run wrapper (assuming reverse_dosimetry() is available in your environment)
#' out <- reverse_dosimetry_wrapping(
#'   Enrichment_data_BMD = df,
#'   duration = 40,
#'   BW = 70,
#'   PFAS_mw = PFAS_mw,
#'   pod = "BMD",
#'   PFAS = c("PFOA","PFOS"),
#'   time = c(1,4),
#'   compartment = c("liver","serum")
#' )
#'
#' # If you add a return() as suggested in @return:
#' # out$reverse_dosimetry_list[["PFOA_1_liver"]]
#' # out$reverse_dosimetry_list_summarized[["PFOS_4_serum"]]
#' }
#'
#' @seealso \code{\link{reverse_dosimetry}}
#'
#' @import PfasDosim
#' @export
#'
reverse_dosimetry_wrapping = function(Enrichment_data_BMD,
                                      duration = 40,
                                      BW = 70,
                                      PFAS_mw = c("PFOA"=414.07,"PFOS"=500.13),
                                      pod = "BMD",
                                      PFAS = c("PFOA","PFOS"),
                                      time = c(1,4,10,14),
                                      compartments = c("liver","serum"),
                                      isStochastic = TRUE,
                                      isParallel = TRUE, 
                                      Nsamples = 10,
                                      n_cores = 5,
                                      free = FALSE
                                      ){
  reverse_dosimetry_list = list()
  reverse_dosimetry_list_summarized = list()
  
  for(pfas in PFAS ){
    for(time_point in time){
      for(compartment in compartments){
        df = Enrichment_data_BMD[Enrichment_data_BMD$PFAS == pfas & Enrichment_data_BMD$time==time_point,]
        df = unique(df[,c("TermID","Ke_description","PFAS",pod)]) # add time and compartment
        
        df = cbind(df,"time" = time_point, "compartment" = compartment)
        colnames(df)[which(colnames(df)==pod)] = "POD"
        colnames(df)[which(colnames(df)=="Ke_description")] = "key_event"
        df$POD = as.numeric(df$POD)
        df = unique(df)
        df = df[,-1]
        ddpp = df[,c("time","POD","compartment","key_event")]
        
        #original POD from omics data multipled by PFAS molecular weigt
        ddpp$POD = ddpp$POD * PFAS_mw[[pfas]]
        
        rd = reverse_dosimetry(chemical = pfas, 
                               BW = BW, 
                               duration = duration, 
                               bmd_df = ddpp,
                               isStochastic = isStochastic,
                               isParallel = isParallel, 
                               Nsamples = Nsamples,
                               n_cores = n_cores,
                               rtol = 1e-5, atol = 1e-5,
                               free = free)
        
        # Store HED estimated for each sample of the stochastic reverse dosimetry model
        X = c()
        for(i in 1:length(rd[[pfas]])){
          X = rbind(X,do.call(cbind,rd[[pfas]][[i]][["bma"]]))
        }
        X = as.data.frame(X)
        X$exposure = as.numeric(X$exposure)
        X$rel_error = as.numeric(X$rel_error)
        X$POD = as.numeric(X$POD)
        
        
        if(isStochastic){
          # Store HED as the average of each sample of the stochastic reverse dosimetry model
          Y = c()
          for(i in 1:length(rd[[pfas]])){
            
            di = do.call(cbind,rd[[pfas]][[i]][["bma"]])
            di = as.data.frame(di)
            di$exposure = as.numeric(di$exposure)
            di$rel_error = as.numeric(di$rel_error)
            di$POD = as.numeric(di$POD)
            
            di$sd = sd(di$exposure)
            di$var = var(di$exposure)
            di$exposure = median(di$exposure)
            di$rel_error = median(di$rel_error)
            di$POD = median(di$POD)
            
            Y = rbind(Y, unique(di))
          }
          Y = as.data.frame(Y)
          # X$exposure = as.numeric(X$exposure)
          # X$rel_error = as.numeric(X$rel_error)
          # X$POD = as.numeric(X$POD)
          reverse_dosimetry_list_summarized[[paste(pfas, time_point, compartment, sep = "_")]] = Y
          
        }
       
        reverse_dosimetry_list[[paste(pfas, time_point, compartment, sep = "_")]] = X
      }
    }
  }
  
  return(list("reverse_dosimetry_list"=reverse_dosimetry_list,
              "reverse_dosimetry_list_summarized"=reverse_dosimetry_list_summarized))
}

