#' Read an excel file
#'
#' Function is a wrapper for readxl::read_excel function.
#' @param filepath A path to the excel file.
#' @param sheetname An integer value to specify which sheet will be read.
#'
#' @return data.table with the contents of the loaded excel file.
#' @export
read_excel <- function(filepath, sheetname = 2){
  snv_file <- readxl::read_excel(filepath, sheet = sheetname);
  return(snv_file);
}


#' Title
#'
#' @param dataset The input dataset containing CNV information
#'
#' @return Modified dataset
#' @export
#'
preprocess_raw_CNV <- function(dataset){
  dataset <- tidyr::separate(dataset, col = Chromosome.Region, into = c("chromosome", "start_position", "end_position"), sep = "[:-]")
  dataset$start_position <- as.numeric(gsub(",", "", dataset$start_position))
  dataset$end_position <- as.numeric(gsub(",", "", dataset$end_position))
  dataset$Event <- factor(dataset$Event, levels = c("Homozygous Copy Loss", "CN Loss", "Allelic Imbalance", "CN Gain", "High Copy Gain"))
  dataset$Event.code <- as.numeric(dataset$Event)
  dataset$Event.code <- dataset$Event.code - 1
  dataset <- dataset |>
    dplyr::select(chromosome, start_position, end_position, Event.code) |>
    dplyr::group_by(chromosome) |>
    dplyr::arrange(start_position, .by_group = TRUE)
  return(dataset)
}

#' Fill missing CNV segments
#'
#' @param cnv_file a source file containing CNV information.
#'
#' @return a modified CNV file
#' @export
#'
fill_missing_cnv_segments <- function(cnv_file){

  chromosome_ends <- data.frame(chromosome = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8','chr9',
                                               'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15','chr16', 'chr17',
                                               'chr18','chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY'),
                                start_position = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067,
                                                   159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
                                                   115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
                                                   59128983, 63025520, 48129895, 51304566, 155270560, 59373566),
                                end_position =  c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067,
                                                  159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
                                                  115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
                                                  59128983, 63025520, 48129895, 51304566, 155270560, 59373566),
                                Event.code = rep(-1, times = 24))

  chromosome_ends_merge <- chromosome_ends[chromosome_ends$chromosome %in% cnv_file$chromosome, ]

  cnv_file_no_overlap <- rbind(cnv_file, chromosome_ends_merge) |>
    dplyr::group_by(chromosome, start_position) |>
    dplyr::filter(!(Event.code == 2 & dplyr::n() > 1)) |>
    dplyr::filter(Event.code!=2)

  temp <- cnv_file_no_overlap |>
    dplyr::group_by(chromosome) |>
    dplyr::arrange(start_position, .by_group = TRUE) |>
    dplyr::mutate(prev_end = dplyr::lag(end_position)) |>
    dplyr::mutate(prev_end = ifelse(is.na(prev_end), 0, prev_end))

  temp2 <- data.frame(chromosome = character(), start_position = numeric(), end_position = numeric(), Event.code = numeric())

  for (i in 1:(nrow(temp))) {
    df <- data.frame(
      chromosome = temp$chromosome[i],
      start_position = temp$prev_end[i],
      end_position = temp$start_position[i],
      Event.code = 2)
    temp2 <- rbind(temp2, df)
  }

  temp2 <- temp2 |>
    dplyr::filter(start_position != end_position)

  final_file <- rbind(temp2, cnv_file_no_overlap) |>
    dplyr::filter(Event.code != -1) |>
    dplyr::group_by(chromosome) |>
    dplyr::arrange(start_position, .by_group = TRUE)

  final_file <- as.data.frame(lapply(final_file, unlist))

  return(final_file)
}
