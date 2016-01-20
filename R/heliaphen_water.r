# Soil functions and variables

# read heliaphen formatted files
#' @export read_heliaphen
read_heliaphen <- function(file, experiment, position, header) {
  
  # read and subset columns
  r <- read.table(file=file, sep="\t", skip=1, stringsAsFactors=FALSE)[,position]
  # add header names
  names(r) <- header
  # subset rows for single experiment
  r <- r %>% filter(grepl(experiment, Plant_Code))
  return(r)
}

