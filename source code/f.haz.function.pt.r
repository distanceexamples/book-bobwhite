#' Evaluating the probability density of observed distances using the hazard-rate model for points
#'
#' The probability density of observed distances for a given distance x1, given scale and 
#' shape parameter sigma1 and shape1 and given truncation distance w is obtained using: 
#' f.haz.function.pt(x1,sigma1,shape1)/integrate(f.haz.function.pt,0,w,sigma1,shape1)$value
#' 
#' @param dis exact radial distance measurement
#' @param sigma scale parameter 
#' @param shape shape parameter

f.haz.function.pt<-function(dis,sigma,shape) {
  f <- dis*(1-exp(-(dis/sigma)^(-shape)))
  f
}
