#' Variance of baseline density
#' 
#' This function calculates the variance of the baseline density estimate (where factor covariates are set
#' to their baseline levels and continuous covariates are set to zero) using the ML parameter estimates
#' and the hessian matrix
#' 
#' @param par maximum likelihood estimates obtained from optimizing the likelihood function
#' @param hessian Hessian matrix obtained from the same optimization
#' @param delta small proportion for calculating small numerical changes in the parameters 

var.Dblml<-function(par,hessian,delta){
int<-par[3]
log.resd<-par[9]
deriv<-matrix(0,1,9)
deriv[1,3]<-(exp(int+(delta*int))*exp(0.5*exp(log.resd)^2)-exp(int-(delta*int))*exp(0.5*exp(log.resd)^2))/(2*delta*int)
deriv[1,9]<-(exp(int)*exp(0.5*exp(log.resd+(delta*log.resd))^2)-exp(int)*exp(0.5*exp(log.resd-(delta*log.resd))^2))/(2*delta*log.resd)
varian<-deriv%*%solve(hessian)%*%t(deriv)
varian
}
