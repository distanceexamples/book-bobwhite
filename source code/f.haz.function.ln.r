#' Evaluating the probability density of observed distances using the hazard-rate model for line transects
#'
#' The probability density of observed distances for a given distance x1, given scale and 
#' shape parameter sigma1 and shape1 and given truncation distance w is obtained using: 
#' f.haz.function.ln(x1,sigma1,shape1)/integrate(f.haz.function.ln,0,w,sigma1,shape1)$value
#' 
#' @param dis exact radial distance measurement
#' @param sigma scale parameter 
#' @param shape shape parameter

f.haz.function.ln<-function(dis,sigma,shape) {
  f <- (1-exp(-(dis/sigma)^(-shape)))
  f
}
