# Plant functions and variables

# read and tidy raw leaf area data
#' @export area_raw
area_raw <- function(experiment, index) {
  # read heliaphen input file
  path <- paste0("data/",experiment,"/phenotype")
  data <- read_excel(paste0(path,"/",experiment,"_phenotype.xlsx"), sheet="area") 
  
  # remove lines corresponding to missing measurements and format dataset
  data <- data %>%
    filter(!(is.na(length) & is.na(width) & is.na(senescence))) %>%
    mutate(
      time=ymd(as.Date(time, origin="1899-12-30"), tz="Europe/Paris"),
      senescence=ifelse(is.na(senescence), 0, senescence)
    )
  
  # compute leaf area from length and width variables
  data <- data %>%
    mutate(area=leaf_size(length, width)) %>% 
    left_join(index)
  
  return(data)
}

# compute individual leaf area with linear interpolation between nodes
#' @export area_nodes
area_nodes <- function(data) {
  
  # filter negative growth rate while keeping senescence information
  # compute total leaf area by discarding senescence : null leaf areas caused by senescence are replaced by last measured values.
  # TODO : replace horrible test by last(area[area > 0], default=0) (bug in dplyr 0.5)
  data_filtered <- data %>%
    group_by(plant_code, leaf) %>%
    mutate(
      area=ifelse(senescence == 0, rate_positive(area), 0),
      area_total=ifelse(
        senescence > 0,
        if (length(tail(area[area > 0], n=1))==0) 0 else tail(area[area > 0], n=1),
        area
      )) 

  # interpolate individual leaf area over stem nodes
  data_interpolated <- data_filtered %>%
    group_by(plant_code, time) %>% 
    do(interpolate_area_node(.)) %>%
    ungroup()
  
  # return interpolated dataset
  return(data_interpolated)
  
}
  

# compute leaf area dynamics over missing timesteps 
#' @export area_dynamics
area_dynamics <- function(data) {
  
  # fit models describing leaf area dynamics 
  data_model <- data %>%
    mutate(day=yday(time)) %>% 
    group_by(plant_code, leaf) %>% 
    nest() %>% 
    mutate(
      model=map(data, model_growth),
      interpolation=map2(data, model, add_interpolation)
    )
  
  # get interpolated data
  data_area <- data_model %>%
    unnest(interpolation) 
  
  return(data_area)
  
}


# logistic model wrapper
#' @export model_logistic
model_logistic <- function(data) {failwith(NA, nls, quiet=TRUE)(area ~ SSlogis(day, Asym, xmid, scal), data=data)}

# linar model wrapper
#' @export model_linear
model_linear <- function(data) {failwith(NA, lm, quiet=TRUE)(area ~ day, data=data)}

# polynomial model wrapper
#' @export model_polynomial
model_polynomial <- function(data) {failwith(NA, lm, quiet=TRUE)(area ~ poly(day, 2), data=data)}

# splines model wrapper
#' @export model_splines
model_splines <- function(data) {failwith(NA, lm, quiet=TRUE)(area ~ splines::bs(day, 3), data=data)}

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
  data %>%
    complete(time=full_seq(time, 24*60*60)) %>% 
    mutate(
      day=yday(time),
      senescence=zoo::na.locf(senescence)
    ) %>% 
    add_predictions(model, var="area") %>% 
    select(time, day, senescence, area) 
}

# filter negative leaf growth rate, possible to identify peaks with cummax(area)-lag(cummax(area)
#' @export rate_positive
rate_positive <- function (x) {
  for (i in 1:length(x)) {
    delta <- x-lag(x)
    x[i] <- if (delta[i] > 0 | is.na(delta[i])) x[i] else x[i-1]
  }
  return(x)
}

# interpolate green and total leaf area
#' @export interpolate_area_node
interpolate_area_node <- function(data){
  
  # test for plants that were not measured
  if(length(data$length[!is.na(data$length)]) == 0) {
    return(NULL)
  } else {
    
    # create senescence coding for specific case (measured value and senescence)
    data <- data %>% mutate(senescence=ifelse(senescence==2 & area > 0, 3, senescence))
    
    # get actual leaf number 
    n <- max(data$leaf)
    
    # linear interpolation
    # TODO code new interpolation functions over nodes (splinefun or [@Keating1992])
    area_total <- with(data, approxfun(leaf, area_total))
    # decode assuming the following encoding for leaf senescence:
    # 0 = no senescence for leaf n and n+1, 1 = leaf n senescent, 2 = both leaf n and n+1 senescent, 3 = only leaf n+1 senescent
    senescence <- function(x) {switch(x+1, c(0,0), c(1,0), c(1,1), c(0,1))}
    
    data_interpolated <- data_frame(
      leaf=1:n,
      senescence=as.vector(mapply(senescence, data$senescence))[1:n],
      area=area_total(1:n)
    ) %>%
      mutate(senescence=ifelse(is.na(senescence), 0, senescence))
    
    # return
    return(data_interpolated)
  } 
}

