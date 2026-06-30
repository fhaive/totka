#' Launch the Shiny Application
#'
#' @param ... arguments to pass to shinyApp.
#' @inheritParams shiny::shinyApp
#'
#' @importFrom shiny shinyApp
#' @export
launch_app <- function(
  onStart = NULL,
  options = list(),
  enableBookmarking = NULL,
  uiPattern = "/",
  ...
) {
  suppressPackageStartupMessages(library(InteractiveComplexHeatmap))
  shinyApp(
    ui = app_ui,
    server = app_server,
    onStart = onStart,
    options = options,
    enableBookmarking = enableBookmarking,
    uiPattern = uiPattern
  )
}
