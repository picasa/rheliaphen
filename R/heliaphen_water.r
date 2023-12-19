# Soil functions and variables

safe_max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=TRUE), NA)

# read heliaphen formatted files
#' @export read_heliaphen
read_heliaphen <- function(file, experiment, position, header) {
  
  # read and subset columns
  r <- utils::read.table(file=file, sep="\t", skip=1, stringsAsFactors=FALSE)[,position]
  # add header names
  names(r) <- header
  # subset rows for single experiment
  r <- r |> dplyr::filter(grepl(experiment, plant_code))
  return(r)
}

# interpolate soil weight, irrigation and FTSW
#' @export interpolate_water_stress
interpolate_water_stress <- function(data, time) {
  # linear interpolations
  w <- with(data, approxfun(time, weight_soil_t, rule = 2:1))
  f <- with(data, approxfun(time, FTSW, rule = 2:1))
  i <- with(data, approxfun(time, irrigation, rule = 2:1))
  
  # return
  output <- data.frame(
    time=time,
    weight_soil_t=w(time),
    FTSW=f(time),
    irrigation=i(time)
  )
  return(output)
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
  list_files <- tidyr::tibble(file = list.files(path=path, full.names=TRUE, pattern="*.csv"))
  
  # read all csv files in target directory, remove duplicate rows, format date and remove row with date out of range
  data_weight_raw <- list_files |>
    dplyr::group_by(file) %>%
    dplyr::do(read_heliaphen(
      file=.$file, experiment,
      position=list_header_position, header=list_header_labels)) |>
    dplyr::ungroup() |>
    dplyr::select(-file) |>
    dplyr::distinct() |> 
    dplyr::mutate(time = lubridate::dmy_hms(time, tz="Europe/Paris")) |>
    dplyr::filter(lubridate::year(time) == lubridate::year(date_start)) 

  # set failed measurements to NA and non irrigated weight to measured weight (instead of 0)
  data_weight_raw <- data_weight_raw |>
    dplyr::mutate(
      weight_t=ifelse(weight_test==0, NA, weight_t),
      weight_irrigated=ifelse(weight_irrigated==0, weight_t, weight_irrigated)
    )
  
  # get pre and post irrigation weight for (normally) irrigated plant 
  data_weight_pre <- data_weight_raw |>
    dplyr::select(plant_code, time, weight=weight_t)
  
  data_weight_post <- data_weight_raw |>
    dplyr::mutate(time=time + lubridate::minutes(1)) |>
    dplyr::select(plant_code, time, weight=weight_irrigated)
 
  # compute added water as (weight_irrigated - weight_t), 0 if failed measurements
  data_irrigation <- data_weight_raw |>
    dplyr::mutate(weight_water = ifelse(is.na(weight_t), 0, weight_irrigated - weight_t) |>
             as.integer()) |> 
    dplyr::select(plant_code, time, weight_water) 
  
  # bind pre and post irrigation weight, remove duplicated measure if no water was added, compute cumulated irrigation  
  data_weight_complete <- dplyr::bind_rows(data_weight_pre, data_weight_post) |>
    dplyr::distinct(plant_code, weight, .keep_all=TRUE) |>
    dplyr::left_join(data_irrigation) |>
    tidyr::replace_na(list(weight_water=0)) |> dplyr::arrange(plant_code, time) |>
    dplyr::group_by(plant_code) |>
    dplyr::mutate(irrigation = cumsum(weight_water))
    
  # add initial weight : maximum weight at starting date, with NA replaced with mean values
  data_weight_init <- data_weight_complete |>
    dplyr::filter(lubridate::day(time) == lubridate::day(date_start)) |>
    dplyr::group_by(plant_code) |>
    dplyr::summarise(weight_0 = safe_max(weight)) |> dplyr::ungroup() |>
    dplyr::mutate(
      weight_0 = replace(
        weight_0, is.na(weight_0), as.integer(mean(weight_0, na.rm=TRUE)))
      )
      
  # add initial weight and index
  data_weight <- data_weight_complete |>
    dplyr::left_join(data_weight_init) |>
    dplyr::select(plant_code, time, weight, weight_0, weight_water, irrigation) |>
    dplyr::left_join(index) |> dplyr::ungroup()
  
  return(data_weight)
  
}


# compute soil water deficit dynamics by plant code
#' @export soil_water_deficit
soil_water_deficit <- function(data, date_start, date_end, weight_dead, weight_hat=44, awc=0.61, timing="daily") {
  
  # FTSW : Fraction of transpirable soil water
  # TTSW : Total transpirable soil water (g)
  # ATSW : Actual transpirable soil water (g)
  
  # TODO: get mean TTSW from transpiration data 
  # TODO: estimate plant weight dynamics and correct pot weight accordingly
  
  list_names <- c(
    "plant_code","time","weight","weight_soil_t","weight_soil_0","FTSW","irrigation",
    "position","line", "column","genotype","treatment","rep"
  )
  
  # correct pot weight according to experiment (pot weight) and design (hats on stressed plants)
  data_weight_soil <- data |>
    dplyr::mutate(
      weight_soil_t = ifelse(
        treatment=="stress", weight-weight_dead-weight_hat, weight-weight_dead),
      weight_soil_0 = ifelse(
        treatment=="stress", weight_0-weight_dead-weight_hat, weight_0-weight_dead)
    )

  # compute FTSW according to fixed TTSW
  data_water_measure <- data_weight_soil |>
    dplyr::mutate(
      TTSW=weight_soil_0 * awc,
      ATSW=weight_soil_t - (weight_soil_0 * (1 - awc))  ,
      FTSW=ATSW/TTSW
    )
  
  switch(
    timing,
    
    # return water deficit computed at timing of measurements
    measure = {
      data_measure <- data_water_measure |> dplyr::select(dplyr::one_of(list_names))
      return(data_measure)
    },
    
    # compute daily pot weight and FTSW by linear interpolation to get regular 24h timesteps
    # water loss is computed :
    # for stressed plants, water loss is the difference in pot weight
    # for control plants, water loss is equal to the irrigation (or difference in pot weight if irrigation is null)
    daily = {
      
      data_daily <- data_water_measure |> 
        dplyr::group_by(plant_code, position, line, column, genotype, treatment, rep) %>%
        dplyr::do(
          purrr::possibly(interpolate_water_stress, data.frame(NULL))
          (., time = seq(date_start, date_end, by ='days'))) |> 
        dplyr::mutate(
          lag_irrigation = irrigation - dplyr::lag(irrigation),
          lag_weight = - (weight_soil_t - dplyr::lag(weight_soil_t))
        ) |> 
        dplyr::mutate(
          water_loss = dplyr::case_when(
            treatment == "control" & lag_weight > lag_irrigation ~ lag_weight,
            treatment == "control" ~ lag_irrigation,
            treatment == "stress" ~ lag_weight,
            TRUE ~ lag_weight
          )) |> 
        dplyr::select(plant_code, time, weight_soil_t, water_loss, FTSW, irrigation, position:rep) |> 
        dplyr::ungroup()
      
      return(data_daily)
    }
  )
}

# get list of harvestable plants during date range, excluding those from previous dates
#' @export plant_harvest
plant_harvest <- function(data, date_start, date_end, threshold=0.1) {
  
  # get list of harvestable stressed plants during date range
  list_harvest <- data |>
    dplyr::filter(treatment=="stress", time >= date_start, time <= date_end, FTSW < threshold) |>
    dplyr::distinct(plant_code, .keep_all=TRUE)
    
  # get plant codes already harvested
  list_harvest_done <- data |>
    dplyr::filter(treatment=="stress", time < date_start, FTSW < threshold) |>
    dplyr::distinct(plant_code, .keep_all=TRUE)
    
  # get common subset of stressed plants
  list_harvest_stress <- dplyr::anti_join(list_harvest, list_harvest_done, by="plant_code")
  
  # get corresponding control plants
  # TODO fuzzy join with time to select appropriate control plant
  list_harvest_control <- dplyr::semi_join(
    data |> 
      dplyr::filter(treatment=="control", time >= date_start, time <= date_end) |> 
      dplyr::distinct(plant_code, .keep_all=TRUE),
    list_harvest_stress,
    by=c("genotype","rep")
  )
  
  table_harvest <- dplyr::bind_rows(list_harvest_stress, list_harvest_control) |> 
    dplyr::arrange(plant_code)
  
  return(table_harvest)
  
}


# write file for field notations from table of harvestable 
#' @export write_heliaphen
write_heliaphen <- function(data, leaf_position = seq(1,35, by=2)) {
  
  # create files for field observations
  file_area <- tidyr::crossing(
    plant_code=data$plant_code,
    leaf=leaf_position
    ) |>
    dplyr::mutate(length=NA, width=NA, senescence=NA)
  
  file_architecture <- dplyr::tibble(
    plant_code=data$plant_code,
    plant_heigth=NA, stem_diameter=NA,
    leaf_number=NA, plant_stage=NA 
  )
  
  return(list(area=file_area, architecture=file_architecture))
  
}


