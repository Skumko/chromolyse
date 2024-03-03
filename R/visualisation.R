#' Create a Circos visualisation
#'
#' @param sv_data somatic variants (SV) source data file. Format of the data should be:
#' @param sv_focus the subset of translocations to focus on: "itx", "ctx" or "all".
#' @param chromosome_selection an optional vector containing a selection of chromosomes to be included in the visualisation.
#' @param cnv_data An optional data source containing CNV information. Used to create a separate CNV track.
#' @param cnv_track_style the type of CNV track visualisation: "full" is default and will display the CNV baseline event (Allelic imbalance, shown in black). Use "trim" to not include this baseline.
#'
#' @export
#'
visualise_circos <- function(sv_data, cnv_data = NULL, chromosome_selection = NULL, sv_focus="all", cnv_track_style="full"){

  #create a color palette used for cluster assignment
  distinctColorPalette <- function(n) {
    rainbow(n, s = 1, v = 1, start = 0, end = 1)
  }
  n <- max(sv_data$Cluster.x, sv_data$Cluster.y)
  color_vector <- stats::setNames(distinctColorPalette(n), c(1:n))

  #define colors for the CNV track
  if (cnv_track_style == "full"){
    cnv_colors=c("blue", "blue","black", "red", "red")

  }
  else if (cnv_track_style == "trim"){
    cnv_colors=c("blue", "blue", "red", "red")
  }
  else{
    stop("Invalid cnv_track_style option")
  }

  #filter SV data based on focus parameter
  if (sv_focus == "itx") {
    soma <- dplyr::filter(sv_data, Chr.x == Chr.y)

  }
  else if(sv_focus == "ctx"){
    soma <- dplyr::filter(sv_data, Chr.x != Chr.y)

  }
  else if (sv_focus == "all"){
    soma <- sv_data
  }
  else{
    stop("Incorrect SV focus choice")
  }

  #filter SV data based on chromosome selection
  if (!is.null(chromosome_selection)) {
    soma <- dplyr::filter(soma, Chr.x %in% chromosome_selection & Chr.y %in% chromosome_selection)
  }


  # create Circos
  if(!is.null(chromosome_selection))
    circlize::circos.initializeWithIdeogram(chromosome.index = paste0("chr", chromosome_selection))
  else
    circlize::circos.initializeWithIdeogram()
  # if specified, use CNV data for its own track
  if(!is.null(cnv_data)){

    if(!is.null(chromosome_selection)){
      cnv_data <- dplyr::filter(cnv_data, chromosome %in% paste0("chr", chromosome_selection))
    }

    if(cnv_track_style == "trim"){
      cnv_data <- dplyr::filter(cnv_data, Event.code != 2)
    }
    cnv_data <- as.data.frame(lapply(cnv_data, unlist))
    cnv_list <- split(cnv_data, cnv_data$Event.code)

    circlize::circos.genomicTrack(cnv_list,track.height = 0.05, bg.lty = 0, stack = TRUE, panel.fun = function(region, value, ...) {
      i = circlize::getI(...)
      circlize::circos.genomicLines(region, value, col = cnv_colors[i], lwd = 1.5,...)
    })
  }

  # add links to represent translocations
  for (i in 1:nrow(soma)) {
    circlize::circos.link(
      paste("chr", soma$Chr.x[i], sep = ""),
      soma$Position.x[i],
      paste("chr", soma$Chr.y[i], sep = ""),
      soma$Position.y[i],
      col = color_vector[soma$Cluster.x[i]]
    )
  }

}
