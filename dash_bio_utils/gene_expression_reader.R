#  Helper function to simplify data import for Clustergram SOFT Files

importSOFT <- function(filepath) {
  geo_data = getGEO(filename = filepath)
  
  geo_table <- Table(geo_data)
  
  row.names(geo_table) <- geo_table$ID_REF
  
  geo_table[1] <- NULL
  
  return(geo_table)
}


