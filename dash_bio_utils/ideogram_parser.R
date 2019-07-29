# 
#  A simple tab-delimited parser for the ideogram, that parses data
#  from the NCBI Genome Ideogram data bank. The function below returns
#  an array, containing the rows of the dataset as strings.
#  NCBI Genome Ideogram data bank: ftp://ftp.ncbi.nlm.nih.gov/pub/gdp/
#  (grab data from here) NCBI Genome Decoration Page:
#  https://www.ncbi.nlm.nih.gov/genome/tools/gdp
#  
#  The parser will convert NCBI Genome Ideogram data to a python list
#  for use with Ideogram.js.
# 
# :param file_location: The location of the file you want to parse, using a relative path.
# :param header_rows: The header rows you want to remove from your dataset.
# :returns: A list containing the NCBI Genome Ideogram data, where each index
# contains a row of the data set as a string.



ncbi_gdp_to_list <- function(file_location="", header_rows = 1) {
  file = readLines(file_location)
  
  cleaned_file <- gsub("\t", " ", file)
  
  dataset_container <- as.list(cleaned_file)

  dataset_container <- dataset_container[header_rows + 1:length(dataset_container) - 1]
  
  return(dataset_container)
}


