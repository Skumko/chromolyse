#' Read an excel file
#'
#' Function is a wrapper for readxl::read_excel function (see\code{readxl::\link[readxl:read_excel]{read_excel function}} for details).
#' @param filepath A path to the excel file.
#' @param sheetname A value to specify which sheet will be read.
#'
#' @return table with the contents of the excel file.
#' @export
readExcel <- function(filepath, sheetname = 2){
  snv_file <- readxl::read_excel(filepath, sheet = sheetname);
  return(snv_file);
}


#' Transform the Nexus standard format file into the required format.
#'
#' The Nexus file format should contain these columns:
#' - `Chromosome.Region`
#' - `Event``
#' - `Length`
#' - `Cytoband`
#' - `X..of.CNV.Overlap`
#' - `Probe.Median`
#' - `X..Heterozygous`
#' - `Probes`
#' - `Count.of.Gene.Symbols`
#'
#' This function transforms the dataset by only selecting columns needed for
#' CNV purposes (mainly visualisation) and using appropriate names:
#' - `chromosome`, specifies which chromosome is concerned,
#' - `start_position`, the starting position of an affected region,
#' - `end_position`, the end position of an affected region,
#' - `event`, a numeric value (0,5) specifying the CNV event.
#'
#' @param dataset The input Nexus format dataset containing CNV information.
#'
#' @return Modified dataset.
#' @export
#'
preprocessRawNexusFile <- function(dataset){
  dataset <- tidyr::separate(dataset, col = Chromosome.Region, into = c("chromosome", "start_position", "end_position"), sep = "[:-]")
  dataset$start_position <- as.numeric(gsub(",", "", dataset$start_position))
  dataset$end_position <- as.numeric(gsub(",", "", dataset$end_position))
  dataset$Event <- factor(dataset$Event, levels = c("Homozygous Copy Loss", "CN Loss", "Allelic Imbalance", "CN Gain", "High Copy Gain"))
  dataset$event <- as.numeric(dataset$Event)
  dataset$event <- dataset$event - 1
  dataset <- dataset |>
    dplyr::select(chromosome, start_position, end_position, event) |>
    dplyr::group_by(chromosome) |>
    dplyr::arrange(start_position, .by_group = TRUE)
  return(dataset)
}

#' Transform the Battenberg standard format file into the required format.
#'
#' The Battenberg file format should contain these columns:
#' - `chr`
#' - `startpos`
#' - `endpos`
#' - `nMaj1_A`
#' - `nMin1_A`
#' - `frac1_A`
#' - `nMaj2_A`
#' - `nMin2_A`
#' - `frac2_A`
#'
#' This function transforms the dataset by only selecting columns needed for
#' CNV purposes (mainly visualisation) and using appropriate names:
#' - `chromosome`, specifies which chromosome is concerned,
#' - `start_position`, the starting position of an affected region,
#' - `end_position`, the end position of an affected region,
#' - `event`, a numeric value (0,5) specifying the CNV event.
#'
#' @param dataset The input Battenberg format dataset containing CNV information.
#'
#' @return Modified dataset.
#' @export
#'
preprocessRawBattenbergFile <- function(dataset){
  colnames(dataset)[colnames(dataset) == "chr"] <- "chromosome"
  dataset$start_position <- as.numeric(dataset$startpos)
  dataset$end_position <- as.numeric(dataset$endpos)
  dataset$event <- as.numeric(dataset$nMaj1_A)
  dataset$event <- dataset$event - 1
  dataset <- dataset |>
    dplyr::select(chromosome, start_position, end_position, event) |>
    dplyr::group_by(chromosome) |>
    dplyr::arrange(start_position, .by_group = TRUE)
  return(dataset)
}

#' Fill missing CNV segments.
#'
#' For a correct CNV Circos visualisation, the CNV track should have no missing segments - ranges of bases where there is no information.
#' If there are no such missing segments, this function returns the original file with no modifications. Otherwise new segments are generated.
#' @param cnv_file a source file containing the CNV information. The format required is:
#' - `chromosome`, specifies which chromosome is concerned,
#' - `start_position`, the starting position of an affected region,
#' - `end_position`, the end position of an affected region,
#' - `event`, a numeric value (0,5) specifying the CNV event.
#' This format can be achieved by using either manual data preprocessing or using one of the functions (`preprocessRawNexusFile`,`preprocessRawBattenbergFile`).
#'
#' @return a modified CNV file
#' @export
#'
fillMissingSegments <- function(cnv_file){

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
                                event = rep(-1, times = 24))

  chromosome_ends_merge <- chromosome_ends[chromosome_ends$chromosome %in% cnv_file$chromosome, ]

  cnv_file_no_overlap <- rbind(cnv_file, chromosome_ends_merge) |>
    dplyr::group_by(chromosome, start_position) |>
    dplyr::filter(!(event == 2 & dplyr::n() > 1)) |>
    dplyr::filter(event!=2)

  temp <- cnv_file_no_overlap |>
    dplyr::group_by(chromosome) |>
    dplyr::arrange(start_position, .by_group = TRUE) |>
    dplyr::mutate(prev_end = dplyr::lag(end_position)) |>
    dplyr::mutate(prev_end = ifelse(is.na(prev_end), 0, prev_end))

  temp2 <- data.frame(chromosome = character(), start_position = numeric(), end_position = numeric(), event = numeric())

  for (i in 1:(nrow(temp))) {
    df <- data.frame(
      chromosome = temp$chromosome[i],
      start_position = temp$prev_end[i],
      end_position = temp$start_position[i],
      event = 2)
    temp2 <- rbind(temp2, df)
  }

  temp2 <- temp2 |>
    dplyr::filter(start_position != end_position)

  final_file <- rbind(temp2, cnv_file_no_overlap) |>
    dplyr::filter(event != -1) |>
    dplyr::group_by(chromosome) |>
    dplyr::arrange(start_position, .by_group = TRUE) |>
    dplyr::filter(start_position < end_position)

  final_file <- as.data.frame(lapply(final_file, unlist))

  return(final_file)
}


#' Filter a structural variation dataset
#'
#' For the purposes of event identification, only certain structural variants are viable.
#' The function filters a dataset to contain only such variants.
#' @param dataset the dataset to filter
#' @param minSupportedReads the threshold value for minimum number of supporting reads for a translocation. Should be provided in the `DP` column.
#' @param minQuality the threshold value for minimum quality of translocations. Should be provided in the `Qual` column.
#'
#' @return the filtered dataset only containing translocations with the required properties.
#' @export
cleanSVDataset <- function(dataset, minSupportedReads = 3, minQuality = 70){
  cleanDataset <- dataset |> dplyr::filter(Type == "CTX" | Type == "ITX") |> dplyr::distinct()
  tryCatch(
    expr = {
      cleanDataset <- cleanDataset |> dplyr::filter(Qual >= minQuality)
    },
    error = function(e){
      message('Could not filter by quality of reads! Make sure `Qual` column is present')
      print(e)
      return(NULL)
    }
  )

  tryCatch(
    expr = {
      cleanDataset <- cleanDataset |> dplyr::filter(DP >= minSupportedReads)
    },
    error = function(e){
      message('Could not filter by supporting reads! Make sure `DP` column is present')
      print(e)
      return(NULL)
    }
  )
  return(cleanDataset)
}
