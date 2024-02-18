#' Create a Circos visualisation
#'
#' @param sv_data somatic variants (SV) source data file. Format of the data should be:
#' @param sv_focus the subset of translocations to focus on: "itx", "ctx" or "all".
#' @param chromosome_selection an optional vector containing a selection of chromosomes to be included in the visualisation.
#' @param cnv_data An optional data source containing CNV information. Used to create a separate CNV track.
#' @param cnv_track_style the type of CNV track visualisation: "trim" is default and will not display the baseline (shown in black). Use "full" to include a baseline.
#'
#' @export
#'
visualise_circos <- function(sv_data, sv_focus="all", chromosome_selection = NULL, cnv_data = NULL, cnv_track_style="trim"){

  #create a color palette used for cluster assignment
  distinctColorPalette <- function(n) {
    rainbow(n, s = 1, v = 1, start = 0, end = 1)
  }
  n <- max(data$Cluster.x, data$Cluster.y)
  color_vector <- setNames(distinctColorPalette(n), c(1:n))

  # define colors for the CNV track
  cnv_colors=c("blue", "blue", "red", "red")

  #filter data based on focus parameter
  if (sv_focus == "itx") {
    soma <- sv_data %>% filter(Chr.x == Chr.y)

  }
  else if(sv_focus == "ctx"){
    soma <- sv_data %>% filter(Chr.x != Chr.y)

  }
  else if (sv_focus == "all"){
    soma <- sv_data
  }
  else{
    stop("Incorrect SV focus choice")
  }

  #filter data based on chromosome selection
  if (!is.null(chromosome_selection)) {
    soma <- soma %>% filter(Chr.x %in% chromosome_selection & Chr.y %in% chromosome_selection)
  }

  # create Circos
  circos.initializeWithIdeogram()
  # if specified, use CNV data for its own track
  if(!is.null(cnv_data)){
    circos.genomicTrack(nexus_adj_list,track.height = 0.05, bg.lty = 0, stack = TRUE, panel.fun = function(region, value, ...) {
      i = getI(...)
      circos.genomicLines(region, value, col = cnv_colors[i], lwd = 1.5,...)
    })
  }

  # add links to represent translocations
  for (i in 1:nrow(all_trans)) {
    circos.link(
      paste("chr", all_trans$Chr.x[i], sep = ""),
      all_trans$Position.x[i],
      paste("chr", all_trans$Chr.y[i], sep = ""),
      all_trans$Position.y[i],
      col = color_vector[all_trans$Cluster.x[i]]
    )
  }

}
