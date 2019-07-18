library(ape)
library(httr)

UniprotQueryBuilder <- setRefClass(
  "UniprotQueryBuilder",
  fields = list(
    query_fields = "list", 
    query_parameters = "list",
    base_url = "character",
    base_query = "character",
    field_separator = "character",
    parameter_separator = "character"
  ),
  methods = list(
    initialize = function(){
      base_url <<- "https://www.uniprot.org/uniprot/"
      base_query <<- "?query=%s"
      field_separator <<- "+AND+"
      parameter_separator <<- "&"
      query_fields <<- list(
        "accession",
        "active",
        "annotation",
        "author",
        "cdantigen",
        "citation",
        "cluster",
        "count",
        "created",
        "database",
        "ec",
        "evidence",
        "existence",
        "family",
        "fragment",
        "gene",
        "gene_exact",
        "goa",
        "host",
        "id",
        "inn",
        "interactor",
        "keyword",
        "length",
        "lineage",
        "mass",
        "method",
        "mnemonic",
        "modified",
        "name",
        "organelle",
        "organism",
        "plasmid",
        "proteome",
        "proteomecomponent",
        "replaces",
        "reviewed",
        "scope",
        "sequence",
        "sequence_modified",
        "source",
        "strain",
        "taxonomy",
        "tissue",
        "web"
      )
      query_parameters <<- list(
        format = list(
          "html",
          "tab",
          "xls",
          "fasta",
          "gff",
          "txt",
          "xml",
          "rdf",
          "list",
          "rss"
        ),
        columns = list(
          "citation",
          "clusters",
          "comments",
          "domains",
          "domain",
          "ec",
          "id",
          "entry name",
          "existence",
          "families",
          "features",
          "genes",
          "go",
          "go-id",
          "interactor",
          "keywords",
          "last-modified",
          "length",
          "organism",
          "organism-id",
          "pathway",
          "protein names",
          "reviewed",
          "sequence",
          "3d",
          "version",
          "virus hosts"
        ),
        sort = list("score"),
        include = list("yes", "no"),
        compress = list("yes", "no"),
        limit = "int",
        offset = "int"
      )
    },
    validate_query_parameters = function(parameters = NULL){
      if (is.null(parameters)){
        parameters = list()
      }
      validated_parameters <- list()
      for (param in names(parameters)){
        if (param %in% names(query_parameters)){
          if (param %in% list("limit", "offset")){
            if (is.numeric(parameters[[param]])){
              validated_parameters[[param]] <- as.character(parameters[[param]])
            }
          } else if (param == "format"){
            if (parameters[[param]] %in% query_parameters[[param]]){
              validated_parameters[[param]] <- parameters[[param]] 
            }
          } else if (param == "columns"){
            column_entry <- unlist(
              strsplit(
                gsub("+", "", parameters[[param]], fixed = TRUE), ","
              )
            )
            set_a <- unlist(unique(query_parameters[[param]]))
            set_b <- unlist(unique(column_entry))
            validated_items <- intersect(set_a, set_b)
            validated_column_entry <- list()
            for (i in 1:length(column_entry)){
              if (column_entry[[i]] %in% validated_items){
                validated_column_entry[[i]] <- column_entry[[i]]
              }
            }
            validated_parameters[[param]] <- paste(
              unlist(validated_column_entry), collapse = ","
            )
          } else if (param %in% list("include", "compress", "sort")){
            if (parameters[[param]] %in% query_parameters[[param]]){
              validated_parameters[[param]]  <- parameters[[param]]
            }
          }
        }
      }
      validated_parameters
    },
    build_query = function(query, fields = NULL, parameters = NULL){
      if (is.null(fields)){
        fields <- list()
      }
      if (is.null(parameters)){
        parameters <- list()
      }
      URL <- paste0(base_url, sprintf(base_query, query))
      for (fieldname in names(fields)){
        if (fieldname %in% query_fields){
          if (is.character(fields[[fieldname]])){
            URL <- paste0(URL, field_separator, fieldname, ":", fields[[fieldname]])
          } else {
            stop(
              sprintf("The value of the field %s is not a string format", fieldname)
            )
          }
        }
      }
      validated_parameters <- validate_query_parameters(parameters)
      for (param in names(validated_parameters)){
        URL <- paste0(
          URL, parameter_separator, param, "=", validated_parameters[[param]]
        )
      }
      URL
    },
    query_into_dataframe = function(
      query, 
      fields = NULL, 
      parameters = NULL, 
      names = NULL
      ){
      target_url <- build_query(query, fields = fields, parameters = parameters)
      print(target_url)
      col_id <- "columns"
      col_names <- NULL
      if (is.null(names)){
        db <- read.csv(URLencode(target_url), sep = "\t")
      } else {
        db <- ape::read.gff(file = URLencode(target_url))
        colnames(db) <- names
      }
      db
    }
  )
) 

pfam_domain_parser <- function(accession){
  URL <- sprintf("http://pfam.xfam.org/protein/%s/graphic", accession)  
  r <- GET(URL)
  jsonData <- content(r, "parsed")
  toJSON(jsonData)
}

#TESTS:
#test <- UniprotQueryBuilder$new()

#test$query_into_dataframe(
  #query = "DDX3X", 
  #fields = list("revieved" = "yes", "database" = "pfam"), 
  #parameters = list(
    #limit = 1, 
    #columns = "id,entry name,length,genes,organism", 
    #sort = "score", 
    #format = "tab"
  #)
#)

#test$query_into_dataframe(
  #query = "DDX3X", 
  #fields = list(revieved = "yes", "database" = "pfam", accession = "O00571"), 
  #parameters = list(format = "gff"), 
  #names = c("name", "db", "mut", "start", "end", "x1", "x2", "x3", "note")
#")
