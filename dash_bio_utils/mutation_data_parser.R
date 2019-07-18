library(jsonvalidate)
library(jsonlite)
library(httr)
# to make loading gff files easier:
library(ape)
library(plyr)
library(dplyr)

### JSON Schemas: ###

EMPTY_MUT_DATA <- list(
  x = list(),
  y = list(),
  mutationGroups = list(),
  domains = list()
)

PFAM_DOM_SCHEMA <-'{
    "type": "array",
    "items": {
        "type": "object",
        "properties": {
            "region": {
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": {
                        "text": {"type": "string"},
                        "start": {"type": "string"},
                        "end": {"type": "string"},
                    },
                    "required": ["text", "start", "end"]
                }
            },
        },
        "required": ["regions"]
    }
}'

PROT_DOM_SCHEMA <- '{
    "type": "array",
    "items": {
        "type": "object",
        "properties": {
            "name": {"type": "string"},
            "coord": {"type": "string"},
        },
        "required": ["name", "coord"]
    }
}'

ARR_SCHEMA <- '{
    "type": "array",
    "items": {"type": ["string", "number"]}
}'

MUT_DATA_SCHEMA <- sprintf('{
  "type": "object",
  "properties": {
    "x": %s,
    "y": %s,
    "mutationGroups": %s,
    "domains": %s,
  },
  "required": ["x"]}', 
  ARR_SCHEMA, ARR_SCHEMA, ARR_SCHEMA, PROT_DOM_SCHEMA
)

# Helper function to convert gff dataframe to a list
# that can be accepted by dashbioNeedlePlot
parse_mutations_uniprot_data <- function(
  gff_data, 
  start = "start", 
  stop = "end", 
  mut_types_to_skip = NULL
  ){
  if (is.null(mut_types_to_skip)){
    mut_types_to_skip <- list(
      "Chain",
      "Region"
    )
  }
  if (!"Chain" %in% mut_types_to_skip){
    mut_types_to_skip <- c(mut_types_to_skip, list("Chain"))
  }
  # Selects the various mutations types in the dataset, except types contained in the above list
  mut_types <- as.character(
    unique(gff_data[!gff_data$mut %in% mut_types_to_skip, "mut"])
  )
  x <- list()
  y <- list()
  mutationgroups <- list()

  for (mut in mut_types){
    data_coord <- gff_data[gff_data$mut == mut, c(start, stop)]
    # split between single and multi-site coordinates
    single_sites <- data_coord[data_coord[,start] == data_coord[,stop],]
    multi_sites <- data_coord[data_coord[,start] != data_coord[,stop],]
    multi_sites[,start] <- paste(
      as.character(multi_sites[,start]),
      as.character(multi_sites[,stop]),
      sep = "-"
    )
    sorted_data <- plyr::count(c(single_sites[,start], multi_sites[,start])) %>% 
      mutate(mut = mut)

    x <- c(x, as.character(sorted_data[,"x"])) 
    y <- c(y, as.character(sorted_data[,"freq"]))
    mutationgroups <- c(mutationgroups, sorted_data[, "mut"])
  }
  # order the results by frequency of occurrence
  order_df <- as.data.frame(
    sort(
      table(unlist(mutationgroups)), 
      decreasing = TRUE
    )
  )
  order_df <- merge(
    data.frame(
      mut = unlist(mutationgroups),
      x = unlist(x),
      y = unlist(y)
    ),
    order_df,
    by.x = "mut", by.y = "Var1" 
  )
  order_df <- order_df[order(order_df[, "Freq"], decreasing = TRUE),]
  
  formatted_data = list(
    "x" = as.character(order_df$x),
    "y" = as.numeric(order_df$y),
    "mutationGroups" = as.character(order_df$mut),
    domains = list()
  )
  formatted_data
}

# helper function to parse protein domains data
# domain_data should adhere to the PFAM_DOM_SCHEMA 
# or the PROT_DOM_SCHEMA defined above
parse_protein_domains_data <- function(domain_data){
  region_key <- "regions"
  region_name_key <- "text"
  region_start_key <- "start"
  region_stop_key <- "end"
  formatted_data <- list()

  if (json_validate(domain_data, PFAM_DOM_SCHEMA)){
    regionlist <- fromJSON(
      domain_data, simplifyVector = FALSE
    )[[1]][[region_key]]
    for (i in 1:length(regionlist)){
      formatted_data[[i]] <- list(
        name = unlist(regionlist[[i]][[region_name_key]]),
        coord = sprintf(
          "%s-%s",
          regionlist[[i]][[region_start_key]],
          regionlist[[i]][[region_stop_key]]
        )
      )
    }
  } else if (json_validate(domain_data, PROT_DOM_SCHEMA)){
    formatted_data <- domain_data
  }
  formatted_data
}

# Helper function to parse mutation data 
parse_mutation_data <- function(mutation_data){
  data <- EMPTY_MUT_DATA 
  mutation_data <- fromJSON(mutation_data)
  for (k in names(data)){
    data[k] <- mutation_data[k]
  }
  data
}

load_protein_domains <- function(accession){
  domain_data <- pfam_domain_parser(accession)
  parse_protein_domains_data(domain_data)
}

load_mutation_data <- function(json_fname = NULL){
  #take a json object and extract the mutation data based one the schema EMPTY_MUT_DATA
  if (!is.null(json_fname)){
    mutationData <- readLines(json_fname)
    return(parse_mutation_data(mutationData))
  } else {
    return(EMPTY_MUT_DATA)
  }
}

parse_mutation_upload_file <- function(contents, fname){
  data <- EMPTY_MUT_DATA
  if (endsWith(fname, ".json"))
    content_string <- unlist(strsplit(contents, ","))[2]
    decoded <- base64_dec(content_string)
    data <- fromJSON(rawToChar(decoded))
  data
}

parse_domain_upload_file <- function(contents, fname){
  data <- list()
  if (endsWith(fname, ".json")){
    content_string <- unlist(strsplit(contents, ","))[2]
    decoded <- base64_dec(content_string)
    data <- fromJSON(rawToChar(decoded))
  }
}

