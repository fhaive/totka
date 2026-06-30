
#' @title Compile Forward Dosimetry Simulation Results
#'
#' @description
#' Takes the nested output of `output_forward_dosimetry()` and compiles,
#' for each PFAS substance, the liver and serum concentration time courses
#' into unified data frames. Each output combines:
#' - model predictions from individual TK models, and  
#' - Bayesian Model Averaging (BMA) posterior draws (renamed with prefixes).
#'
#' The structure of the nested list looks like this: 
#'   
#   output_forward_dosimetry
# ├── PFOA
# │    ├── bma
# │    │     ├── serum_df   (16 × 4001 matrix-like df)
# │    │     └── liver_df   (16 × 4001)
# │    └── single_models
# │          ├── serum_df   (16 × 5)
# │          └── liver_df   (16 × 4)
# └── PFOS
# ├── bma
# │     ├── serum_df   (16 × 4001)
# │     └── liver_df   (16 × 4001)
# └── single_models
# ├── serum_df   (16 × 5)
# └── liver_df   (16 × 4)
#' 
#' @details
#' This function standardizes the structure of dosimetry simulation output by:
#' 1. Extracting time columns and single‑model predictions  
#' 2. Renaming BMA draw columns using a consistent prefix (e.g. `"Cliver_BMA"`)  
#' 3. Combining time + single models + BMA draws into a single dataframe  
#'
#' The output is a list with one entry per PFAS, each containing:
#' - `liver`: concentration‑time profiles in liver  
#' - `serum`: concentration‑time profiles in serum  
#'
#' @param output_forward_dosimetry
#' A named list (usually from `output_forward_dosimetry()`) where each
#' element represents a PFAS (e.g., `"PFOA"`, `"PFOS"`), containing:
#' - `bma$liver_df`, `bma$serum_df`: BMA posterior draws  
#' - `single_models$liver_df`, `single_models$serum_df`: predictions from
#'   individual TK models  
#'
#' @return
#' A named list of PFAS-specific lists with elements:
#' \describe{
#'   \item{liver}{A dataframe with time, single-model predictions, and BMA draws}
#'   \item{serum}{A dataframe with time, single-model predictions, and BMA draws}
#' }
#'
#' @examples
#' \dontrun{
#' out <- output_forward_dosimetry(...)
#' compiled <- compiling_forward_dosimetry_simulation_dataframes_all_samples(out)
#'
#' names(compiled)
#' compiled$PFOA$liver
#' compiled$PFOS$serum
#' }
#'
#' @export
compiling_forward_dosimetry_simulation_dataframes_all_samples <- function(output_forward_dosimetry) {
  result <- list()
  
  for (pfas in names(output_forward_dosimetry)) {
    res <- output_forward_dosimetry[[pfas]]
    
    liver_df <- build_combined_df(
      bma_df    = res$bma$liver_df,
      single_df = res$single_models$liver_df[, -1, drop = FALSE],
      prefix    = "Cliver_BMA"
    )
    
    serum_df <- build_combined_df(
      bma_df    = res$bma$serum_df,
      single_df = res$single_models$serum_df[, -1, drop = FALSE],
      prefix    = "Cserum_BMA"
    )
    
    result[[pfas]] <- list(
      liver = liver_df,
      serum = serum_df
    )
  }
  
  return(result)
}


