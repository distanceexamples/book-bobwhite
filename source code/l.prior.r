#' Uniform prior probabilities
#' 
#' This function returns the summed log of the uniform prior probabilities 
#' @param coef vector of coefficient values
#' @param lo vector of lower boundaries for uniform distributions
#' @param up vector of upper boundaries for uniform distributions 
#' 
#' @details 
#' In contrast to using *sum(log(dunif(...)))*, the \code{l.prior} function returns 
#' a very small value for the log of the prior probability if the value is 
#' outside the boundaries of the distribution. This prevents the acceptance of 
#' a newly proposed parameter value if it is outside the boundaries. 
#' In contrast, using *log(dunif(...))* returns *-Inf* if the parameter value is outside 
#' the boundaries which may cause the algorithm to crash. 
#' 
l.prior<-function(coef,lo,up){
  l.u.int<-log(dunif(coef,lo,up))
  if(any(abs(l.u.int)==Inf))l.u.int<- -100000
  sum(l.u.int)
}
