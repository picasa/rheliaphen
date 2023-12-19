# Plant functions and variables

# read and tidy raw leaf area data
#' @export area_raw
area_raw <- function(experiment, index, method = "manual") {
  
  # read heliaphen input file
  path <- paste0("data/",experiment,"/phenotype")

  switch(
    method,
    
    # leaf measurements include length and width
    manual = {
      
      data <- readxl::read_excel(paste0(path,"/",experiment,"_phenotype.xlsx"), sheet="area") 
      
      # remove lines corresponding to missing measurements and format dataset
      data <- data |>
        dplyr::filter(!(is.na(length) & is.na(width) & is.na(senescence))) |>
        dplyr::mutate(
          time = lubridate::ymd(as.Date(time, origin="1899-12-30"), tz="Europe/Paris"),
          senescence = ifelse(is.na(senescence), 0, senescence)
        )
      
      # compute leaf area from length and width variables
      data <- data |>
        dplyr::mutate(area = leaf_size(length, width)) |> 
        dplyr::left_join(index)
    },
    
    # leaf measurements include only area
    sensor = {
      
      data <- readxl::read_excel(paste0(path,"/",experiment,"_phenotype.xlsx"), sheet="sensor") 
      
      # remove lines corresponding to missing measurements and format dataset
      data <- data |>
        dplyr::filter(!(is.na(area) & is.na(senescence))) |>
        dplyr::mutate(
          time = lubridate::ymd(as.Date(time, origin="1899-12-30"), tz="Europe/Paris"),
          senescence = ifelse(is.na(senescence), 0, senescence)
        )
      
      data <- data |> dplyr::left_join(index)
    }
  )
  
  return(data)
}

# compute individual leaf area with linear interpolation between nodes
#' @export area_nodes
area_nodes <- function(data) {
  
  # filter negative growth rate while keeping senescence information
  # compute total leaf area by discarding senescence : null leaf areas caused by senescence are replaced by last measured values.
  
  data_filtered <- data |>
    dplyr::group_by(plant_code, leaf) |>
    dplyr::mutate(
      area = ifelse(senescence == 0, rate_positive(area), 0),
      area_total = ifelse(
        senescence > 0,
        if (length(utils::tail(area[area > 0], n=1))==0) 0 else utils::tail(area[area > 0], n=1),
        area
      )) 

  # interpolate individual leaf area over stem nodes
  data_interpolated <- data_filtered |>
    dplyr::group_by(plant_code, time) %>%
    dplyr::do(interpolate_area_node(.)) |>
    dplyr::ungroup()
  
  # return interpolated dataset
  return(data_interpolated)
  
}
  

# compute leaf area dynamics over missing timesteps 
#' @export area_dynamics
area_dynamics <- function(data) {
  
  # fit models describing leaf area dynamics 
  data_model <- data |>
    dplyr::mutate(day = lubridate::yday(time)) |> 
    dplyr::group_by(plant_code, leaf) |> 
    tidyr::nest() |> 
    dplyr::mutate(
      model = purrr::map(data, model_growth),
      interpolation = purrr::map2(data, model, add_interpolation)
    )
  
  # get interpolated data
  data_area <- data_model |>
    tidyr::unnest(interpolation, .drop=TRUE) 
  
  return(data_area)
  
}


# predict plant leaf area from image analysis
#' @export area_predict
area_predict <- function(experiment, index, model, min = 0.5) {
  
  # read features from ipsophen
  data_area_image <- readr::read_csv(glue::glue("data/{experiment}/phenotype/{experiment}_ipsophen.csv"))
  
  # shape raw data
  data_area_raw <- data_area_image |> 
    dplyr::mutate(
      time = date_time,
      day = lubridate::yday(time),
      experiment = stringr::str_to_upper(experiment),
      position = stringr::str_to_upper(plant)) |> 
    dplyr::left_join(index, by = c("experiment", "position")) |> 
    dplyr::select(
      plant_code, treatment, genotype, rep, time, day,
      area_sens = area, hull_area, shape_height) 
  
  # predict plant area from selected features using a linear model, average per day
  data_prediction_raw <- data_area_raw |> modelr::add_predictions(model) |> 
    dplyr::filter(pred > 0) |> 
    dplyr::group_by(plant_code, treatment, genotype, rep, time = as.Date(time), day) |> 
    dplyr::summarise(area = mean(pred), area_sd = stats::sd(pred)) 
  
  # post-process predictions : remove plant with low growth (less than 1.5 initial size) and smooth dynamics with spline 
  data_prediction_post <- data_prediction_raw |> 
    dplyr::group_by(plant_code) |> tidyr::nest() |> 
    dplyr::mutate(growth = purrr::map_dbl(data, ~ (max(..1$area) - first(..1$area)) / first(..1$area))) |> 
    dplyr::filter(growth > min) |> 
    dplyr::mutate(
      model = purrr::map(data, model_splines),
      prediction = purrr::map2(data, model, modelr::add_predictions, var="area_splines")
    ) |>
    dplyr::select(-data, -model, -growth) |> tidyr::unnest(prediction) |> 
    dplyr::mutate(
      time = lubridate::ymd(time, tz="Europe/Paris"),
      area = area_splines) |>
    dplyr::select(plant_code:area)
    
  return(data_prediction_post)
  
}


# logistic model wrapper
#' @export model_logistic
model_logistic <- function(data) {
  purrr::possibly(stats::nls, NA)(area ~ SSlogis(day, Asym, xmid, scal), data=data)
}

# linar model wrapper
#' @export model_linear
model_linear <- function(data) {purrr::possibly(stats::lm, NA)(area ~ day, data=data)}

# polynomial model wrapper
#' @export model_polynomial
model_polynomial <- function(data) {purrr::possibly(stats::lm, NA)(area ~ stats::poly(day, 2), data=data)}

# splines model wrapper
#' @export model_splines
model_splines <- function(data) {purrr::possibly(stats::lm, NA)(area ~ splines::bs(day, 3), data=data)}

# growth model wrapper
#' @export model_growth
model_growth <- function(data) {
  # fit for default logictic model
  fit_nonlinear <- model_logistic(data)
  # fit for polynomial or linear model
  if (is.na(fit_nonlinear[1])) {
    fit_linear <- if (nrow(data) > 2) model_polynomial(data) else model_linear(data)
  } 
  # select fit to return
  return(if (is.na(fit_nonlinear[1])) fit_linear else fit_nonlinear)
}

# add interpolated data as a function of data and model
#' @export add_interpolation
add_interpolation <- function(data, model) {
  
  # complete dataframe with missing timesteps 
  data |>
    tidyr::complete(time = tidyr::full_seq(time, 24*60*60)) |> 
    dplyr::mutate(
      day = lubridate::yday(time),
      senescence=zoo::na.locf(senescence)
    ) |> 
    modelr::add_predictions(model, var="area") |> 
    dplyr::select(time, day, senescence, area) 
}

# filter negative leaf growth rate, possible to identify peaks with cummax(area)-lag(cummax(area)
#' @export rate_positive
rate_positive <- function (x) {
  for (i in 1:length(x)) {
    delta <- x - dplyr::lag(x)
    x[i] <- if (delta[i] > 0 | is.na(delta[i])) x[i] else x[i-1]
  }
  return(x)
}

# interpolate green and total leaf area
#' @export interpolate_area_node
interpolate_area_node <- function(data){
  
  # test for plants that were not measured
  if(length(data$length[!is.na(data$length)]) == 0) {
    return(data.frame(NULL))
  } else {
    
    # create senescence coding for specific case (measured value and senescence)
    data <- data |> dplyr::mutate(senescence = ifelse(senescence==2 & area > 0, 3, senescence))
    
    # get actual leaf number 
    n <- max(data$leaf)
    
    # linear interpolation
    # TODO code new interpolation functions over nodes (splinefun or [@Keating1992])
    area_total <- with(data, stats::approxfun(leaf, area_total))
    # decode assuming the following encoding for leaf senescence:
    # 0 = no senescence for leaf n and n+1, 1 = leaf n senescent, 2 = both leaf n and n+1 senescent, 3 = only leaf n+1 senescent
    senescence <- function(x) {switch(x+1, c(0,0), c(1,0), c(1,1), c(0,1))}
    
    data_interpolated <- dplyr::data_frame(
      leaf = 1:n,
      senescence = as.vector(mapply(senescence, data$senescence))[1:n],
      area = area_total(1:n)
    ) |>
      dplyr::mutate(senescence = ifelse(is.na(senescence), 0, senescence))
    
    # return
    return(data_interpolated)
  } 
}

