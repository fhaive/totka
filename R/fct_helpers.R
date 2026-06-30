
#' get_best_match 
#'
#' @description Get the best match for string in the vector
#'
#' @return Returns best string match or NULL.
#'
#' @importFrom dplyr mutate arrange pull
#' @importFrom stringr str_detect regex str_length
#' @importFrom purrr keep pluck
#' @importFrom tibble as_tibble
#' @export
#' @noRd
get_best_match <- function(x, y){
    res <- x |> purrr::keep(\(x){stringr::str_detect(x, pattern=stringr::regex(y, ignore_case=TRUE))})
    if(length(x)==0){
        res <- NULL
    }else{
        res <- res |> tibble::as_tibble() |> dplyr::mutate(len=stringr::str_length(value)) |> dplyr::arrange(len) |> dplyr::pull(value) |> purrr::pluck(1)
    }
    return(res)
}

#' sort_alphnum 
#'
#' @description Sort vector of alphanumeric strings str_num 
#'
#' @return Returns sorted alphanumeric strings vector.
#'
#' @import dplyr
#' @importFrom stringr str_split
#' @importFrom purrr map
#' @importFrom tibble as_tibble
#' @export
#' @noRd
sort_alphnum <- function(x){
  res <- x |> unique() |> stringr::str_split(pattern="_") |> purrr::map(\(x){c(paste0(x[1:length(x)-1],collapse="_"), x[length(x)])}) |> as.data.frame() |> t() |> as.data.frame() |> tibble::as_tibble() |> dplyr::mutate(label=paste0(V1,"_",V2)) |> dplyr::mutate(V2=V2 |> as.numeric()) |> dplyr::arrange(V1,V2) |> dplyr::pull(label)
  return(res)
}

# =====================================================================
# Internal Helper Functions (Roxygen-tagged but not exported)
# =====================================================================

#' @title Extract and Rename BMA Draw Columns
#' @description Internal utility to rename BMA posterior draw columns.
#' @param df Dataframe from BMA output; first column must be time.
#' @param prefix Prefix for renamed draw columns (e.g. `"Cliver_BMA"`).
#' @return A dataframe with renamed posterior draw columns.
#' @keywords internal
extract_bma_draws <- function(df, prefix) {
  draws <- df[, -1, drop = FALSE]
  colnames(draws) <- paste(prefix, colnames(draws), sep = "_")
  draws
}

#' @title Construct Combined Simulation Dataframe
#' @description
#' Internal helper to merge time, single-model predictions,
#' and renamed BMA draws into one dataframe.
#' @param bma_df BMA dataframe (first column = time).
#' @param single_df Dataframe of single-model predictions.
#' @param prefix Prefix for BMA draws.
#' @return A unified simulation dataframe.
#' @keywords internal
build_combined_df <- function(bma_df, single_df, prefix) {
  
  time_col <- bma_df[[1]]
  bma_draws <- extract_bma_draws(bma_df, prefix)
  
  combined <- cbind(
    time = time_col,
    single_df,
    bma_draws
  )
  
  combined
}


