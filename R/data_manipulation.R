#' Read an excel file
#'
#' Function is a wrapper for readxl::read_excel function.
#' @param filepath A path to the excel file.
#' @param sheetname An integer value to specify which sheet will be read.
#'
#' @return data.table with the contents of the loaded excel file.
#' @examples
#' temp1 <- read_excel("./my_excel_file.xlsx");
#' temp1 <- read_excel("./my_excel_file.xls", sheetname=3);
#' @export
read_excel <- function(filepath, sheetname = 2){
  snv_file <- readxl::read_excel(filepath, sheet = sheetname);
  return(snv_file);
}
