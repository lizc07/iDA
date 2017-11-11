#' @include server.R
#' @include ui.R
#' @include global.R
NULL

#'
#' iteractive Data Analysis
#'
#' Interactive Data Analysis for scRNA-Seq using tag-based UMI data
#'
#' @param NULL
#'
#' @export
#'
#' @examples
#' # runApp
#' # iDA()
iDA <- function(){
  options(shiny.maxRequestSize=1000*1024^2)
  pbmc <<- NULL
  pbmc.markers <<- NULL

  pal_heatmap1 <<- colorpanel(100,"blue","white","red")
  pal_heatmap2 <<- colorpanel(100,"darkblue","white","red")

  runApp(shinyApp(ui = ui, server = server))
}
