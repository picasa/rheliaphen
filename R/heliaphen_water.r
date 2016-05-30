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

# compute soil weight dynamics by plant code
# TODO: set columns and default path
#' @export soil_weight
soil_weight <- function(experiment, index, date_start) {
  
  # list column position to keep in csv files
  list_header_position <- c(1,6,7,9,11)
  # list header labels for these columns
  list_header_labels <- c("Plant_Code","Weight_test","Weight_t","Date","Weight_irrigated")
  
  # list files to read
  path <- paste0("data/",experiment,"/raw")
  list_files <- data_frame(file=list.files(path=path, full.names=TRUE, pattern="*.csv"))
  
  # read all csv files in target directory, remove duplicate rows, format date and remove row with date out of range
  data_weight <- list_files %>%
    group_by(file) %>%
    do(read_heliaphen(file=.$file, experiment, position=list_header_position, header=list_header_labels)) %>%
    ungroup() %>%
    select(-file) %>%
    distinct() %>% 
    mutate(Date=dmy_hms(Date)) %>%
    filter(year(Date)==year(date_start)) 

  # set failed measurements to NA
  # TODO check Weight_irrigated==0 : (1) weight after irrigation or (2) weight of added water
  data_weight <- data_weight %>%
    mutate(
      Weight_t=ifelse(Weight_test==0, NA, Weight_t),
      Weight_irrigated=ifelse(Weight_irrigated==0, Weight_t, Weight_irrigated)
    )
  
  # add initial weight, replacing eventual NAs with mean initial value, then use max value of Weight_t for initial weight
  # TODO: check if row number == plant number
  data_weight_init <- data_weight %>%
    filter(day(Date)==day(date_start)) %>%
    mutate(Weight_0=replace(Weight_t, is.na(Weight_t), as.integer(mean(Weight_t, na.rm=TRUE)))) %>%
    group_by(Plant_Code) %>%
    summarise(Weight_0=max(Weight_0)) %>%
    select(Plant_Code, Weight_0) 
  
  # set Weight_irrigated to Weight_0 at initial date
  data_weight <- data_weight %>%
    left_join(data_weight_init) %>%
    mutate(Weight_irrigated=ifelse(day(Date)==day(date_start), Weight_0, Weight_irrigated))
  
  # compute plant irrigation as cumsum of (weight_irrigated - weight_t), zero if failed measurements
  data_weight <- data_weight %>%
    group_by(Plant_Code) %>%
    mutate(
      Weight_water=as.integer(ifelse(is.na(Weight_t), 0, Weight_irrigated - Weight_t)),
      Irrigation=cumsum(Weight_water)
    )
  
  # add index
  data_weight <- data_weight %>% left_join(index) %>% select(-Weight_test)
  
  return(data_weight)
}


# compute soil water deficit dynamics by plant code
#' @export soil_water_deficit
soil_water_deficit <- function(data, index, date_start, date_end, weight_dead=450, awc=0.39, timing="daily") {
  
  list_names <- c(
    "Plant_Code","Date","Weight_t","FTSW","Irrigation",
    "Position","Line", "Column","Genotype","Treatment","Repeat"
  )
  
  # compute FTSW according to fixed TTSW
  # TODO account for equipment of stressed pot
  # TODO: get mean TTSW from transpiration data 
  # TODO: use plant weight at last date as lower TTSW bound)
  # TODO: estimate plant weight and adjust pot weight
  data_ftsw_measure <- data %>%
    mutate(
      ATSW=(Weight_t-weight_dead) - (Weight_0-weight_dead) * awc,
      TTSW=(Weight_0-weight_dead) - (Weight_0-weight_dead) * awc,
      FTSW=ATSW/TTSW
    )
  
  switch(
    timing,
    
    # return water deficit computed at timing of measurements
    measure = {
      data_ftsw_measure <- data_ftsw_measure %>% select(one_of(list_names))
      
      return(data_ftsw_measure)
    },
    
    # compute daily pot weight and FTSW by linear interpolation to get regular 24h timesteps
    daily = {
      
      data_ftsw_daily <- data_ftsw_measure %>%
        ddply(
          ~ Plant_Code,
          plyr::failwith(NULL, interpolate_water_stress),
          time=seq(date_start, date_end, by='days')
        )
      
      # add index
      data_ftsw_daily <- data_ftsw_daily %>%
        left_join(index) %>%
        select(one_of(list_names))
        
      return(data_ftsw_daily)
    }
  )
}

# get list of harvestable plants during date range, excluding those from previous dates
#' @export plant_harvest
plant_harvest <- function(data, date_start, date_end, threshold=0.1) {
  
  # get list of harvestable stressed plants during date range
  list_harvest <- data %>%
    filter(Treatment=="stress", Date >= date_start, Date <= date_end, FTSW < threshold) %>%
    group_by(Plant_Code) %>%
    distinct()
  
  # get plant codes already harvested
  list_harvest_done <- data %>%
    filter(Treatment=="stress", Date < date_start, FTSW < threshold) %>%
    group_by(Plant_Code) %>%
    distinct()
  
  list_harvest_stress <- anti_join(list_harvest, list_harvest_done, by="Plant_Code")
  
  # get corresponding control plants
  list_harvest_control <- semi_join(
    data %>% filter(Treatment=="control"),
    list_harvest_stress,
    by=c("Date","Genotype","Repeat")
  )
  
  table_harvest <- bind_rows(list_harvest_stress, list_harvest_control) %>% arrange(Plant_Code)
  
  return(table_harvest)
  
}


# write file for field notations from table of harvestable 
#' @export write_heliaphen
write_heliaphen <- function(data, leaf_position=seq(1,35, by=2)) {
  
  # create files for field observations
  file_area <- tidyr::crossing(
    Plant_Code=data$Plant_Code,
    leaf=leaf_position
  ) %>%
    mutate(length=NA, width=NA, senescence=NA)
  
  file_architecture <- data_frame(
    Plant_Code=data$Plant_Code,
    plant_heigth=NA, stem_diameter=NA,
    leaf_number=NA, plant_stage=NA 
  )
  
  return(list(area=file_area, architecture=file_architecture))
  
}


