# Soil functions and variables

# read heliaphen formatted files
#' @export read_heliaphen
read_heliaphen <- function(file, experiment, position, header) {
  
  # read and subset columns
  r <- read.table(file=file, sep="\t", skip=1, stringsAsFactors=FALSE)[,position]
  # add header names
  names(r) <- header
  # subset rows for single experiment
  r <- r %>% filter(grepl(experiment, plant_code))
  return(r)
}

# compute soil weight dynamics by plant code
# TODO: set columns and default path
#' @export soil_weight
soil_weight <- function(experiment, index, date_start) {
  
  # list column position to keep in csv files
  list_header_position <- c(1,6,7,9,11)
  # list header labels for these columns
  list_header_labels <- c("plant_code","weight_test","weight_t","time","weight_irrigated")
  
  # list files to read
  path <- paste0("data/",experiment,"/raw")
  list_files <- data_frame(file=list.files(path=path, full.names=TRUE, pattern="*.csv"))
  
  # read all csv files in target directory, remove duplicate rows, format date and remove row with date out of range
  data_weight_raw <- list_files %>%
    group_by(file) %>%
    do(read_heliaphen(file=.$file, experiment, position=list_header_position, header=list_header_labels)) %>%
    ungroup() %>%
    select(-file) %>%
    distinct() %>% 
    mutate(time=dmy_hms(time, tz="Europe/Paris")) %>%
    filter(year(time)==year(date_start)) 

  # set failed measurements to NA and non irrigated weight to measured weight (instead of 0)
  data_weight_raw <- data_weight_raw %>%
    mutate(
      weight_t=ifelse(weight_test==0, NA, weight_t),
      weight_irrigated=ifelse(weight_irrigated==0, weight_t, weight_irrigated)
    )
  
  # get pre and post irrigation weight for (normally) irrigated plant 
  data_weight_pre <- data_weight_raw %>%
    select(plant_code, time, weight=weight_t)
  
  data_weight_post <- data_weight_raw %>%
    mutate(time=time + minutes(1)) %>%
    select(plant_code, time, weight=weight_irrigated)
 
  # compute added water as (weight_irrigated - weight_t), 0 if failed measurements
  data_irrigation <- data_weight_raw %>%
    mutate(weight_water=as.integer(ifelse(is.na(weight_t), 0, weight_irrigated - weight_t))) %>% 
    select(plant_code, time, weight_water) 
  
  # bind pre and post irrigation weight, remove duplicated measure if no water was added, compute cumulated irrigation  
  data_weight_complete <- bind_rows(data_weight_pre, data_weight_post) %>%
    group_by(plant_code, weight) %>% distinct() %>%
    left_join(data_irrigation) %>% replace_na(list(weight_water=0)) %>% arrange(plant_code, time) %>%
    group_by(plant_code) %>%
    mutate(irrigation=cumsum(weight_water))
    
  # add initial weight : maximum weight at starting date, with NA replaced with mean values
  data_weight_init <- data_weight_complete %>%
    filter(day(time)==day(date_start)) %>%
    group_by(plant_code) %>%
    summarise(weight_0=max(weight, na.rm=TRUE)) %>% ungroup() %>%
    mutate(weight_0=replace(weight_0, is.na(weight_0), as.integer(mean(weight_0, na.rm=TRUE))))
      
  # add initial weight and index
  data_weight <- data_weight_complete %>%
    left_join(data_weight_init) %>%
    select(plant_code, time, weight, weight_0, weight_water, irrigation) %>%
    left_join(index) %>% ungroup()
  
  return(data_weight)
  
}


# compute soil water deficit dynamics by plant code
#' @export soil_water_deficit
soil_water_deficit <- function(data, index, date_start, date_end, weight_dead, weight_hat=44, awc=0.39, timing="daily") {
  
  list_names <- c(
    "plant_code","time","weight","FTSW","irrigation",
    "position","line", "column","genotype","treatment","repeat"
  )
  
  # correct pot weight according to experiment and design
  data_weight_soil <- data %>%
    mutate(
      weight_soil_t=ifelse(treatment=="stress", weight-weight_dead-weight_hat, weight-weight_dead),
      weight_soil_0=ifelse(treatment=="stress", weight_0-weight_dead-weight_hat, weight_0-weight_dead)
    )
  
  # compute FTSW according to fixed TTSW
  # TODO account for equipment of stressed pot
  # TODO: get mean TTSW from transpiration data 
  # TODO: use plant weight at last date as lower TTSW bound)
  # TODO: estimate plant weight and adjust pot weight
  data_ftsw_measure <- data_weight_soil %>%
    mutate(
      TTSW=weight_soil_0 - weight_soil_0 * awc,
      ATSW=weight_soil_t - weight_soil_0 * awc,
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
          ~ plant_code,
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
  
  # test format of date
  
  
  # get list of harvestable stressed plants during date range
  list_harvest <- data %>%
    filter(treatment=="stress", time >= date_start, time <= date_end, FTSW < threshold) %>%
    group_by(plant_code) %>%
    distinct()
  
  # get plant codes already harvested
  list_harvest_done <- data %>%
    filter(treatment=="stress", time < date_start, FTSW < threshold) %>%
    group_by(plant_code) %>%
    distinct()
  
  
  # common part
  list_harvest_stress <- anti_join(list_harvest, list_harvest_done, by="plant_code")
  
  # get corresponding control plants
  list_harvest_control <- semi_join(
    data %>% filter(treatment=="control"),
    list_harvest_stress,
    by=c("time","genotype","repeat")
  )
  
  table_harvest <- bind_rows(list_harvest_stress, list_harvest_control) %>% arrange(plant_code)
  
  return(table_harvest)
  
}


# write file for field notations from table of harvestable 
#' @export write_heliaphen
write_heliaphen <- function(data, leaf_position=seq(1,35, by=2)) {
  
  # create files for field observations
  file_area <- tidyr::crossing(
    plant_code=data$plant_code,
    leaf=leaf_position
    ) %>%
    mutate(length=NA, width=NA, senescence=NA)
  
  file_architecture <- data_frame(
    plant_code=data$plant_code,
    plant_heigth=NA, stem_diameter=NA,
    leaf_number=NA, plant_stage=NA 
  )
  
  return(list(area=file_area, architecture=file_architecture))
  
}


