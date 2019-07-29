
#  Helper function to simplify data import for Speck
importSpeck <- function(filepath, 
                        header = FALSE, 
                        skip = 2) {
  textdata <- read.table(
    text = paste0(
      readLines(filepath), collapse="\n"
    ),
    header = header,
    skip = skip,
    col.names = c("symbol", "x", "y", "z"),
    stringsAsFactors = FALSE)
  return(dashTable::df_to_list(textdata))
}