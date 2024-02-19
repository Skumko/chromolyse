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


#' Fill missing CNV segments
#'
#' @param cnv_file a source file containing CNV information.
#'
#' @return a modified CNV file
#' @export
#'
fill_missing_cnv_segments <- function(cnv_file){

}
