#' Evaluating the probability density of observed distances using the half-normal model for point transects
#'
#' The probability density of observed distances for a given distance x1, given scale 
#' parameter sigma1 and given truncation distance w is obtained using: 
#' f.hn.function.pt(x1,sigma1)/integrate(f.hn.function.pt,0,w,sigma1)$value
#' 
#' @param dis exact radial distance measurement
#' @param sigma scale parameter 

f.hn.function.pt<-function(dis,sigma) {
  f <- dis*exp(-dis^2/(2*sigma^2))
  f
}
