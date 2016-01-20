# Plant functions and variables


# interpolate green and total leaf area
#' @export interpolate_leaf_area
interpolate_leaf_area <- function(x){
  
  # test for plants that were not measured
  if(length(x$length[!is.na(x$length)]) == 0) {
    return(NULL)
  } else {
    # get actual leaf number 
    # n <- with(x, leaf[match(last(area[!is.na(area)]), area)]) # index for last non NA area
    n <- max(x$leaf)
    # linear interpolation
    s <- with(x, approxfun(leaf, area))
    g <- with(x, approxfun(leaf, area_total))
    # TODO code senescence interpolation
    
    r <- data.frame(
      leaf=1:n,
      area=s(1:n),
      area_total=g(1:n)
    ) %>% 
    left_join(x %>% select(leaf, senescence))
    
    # TODO : improve syntax
    r$area[which(r$senescence==2)+1] <- 0
    r$area_total[which(r$senescence==2)+1] <- 0
    
    # return
    return(r)
  } 
}

# interpolate weight, irrigation and FTSW
#' @export interpolate_water_stress
interpolate_water_stress <- function(x, time) {
  # linear interpolations
  w <- with(x, approxfun(Date, Weight_t, rule = 2:1))
  f <- with(x, approxfun(Date, FTSW, rule = 2:1))
  i <- with(x, approxfun(Date, Irrigation, rule = 2:1))
  
  # return
  r <- data.frame(
    Date=time,
    Weight_t=w(time),
    FTSW=f(time),
    Irrigation=i(time)
  )
  return(r)
}

