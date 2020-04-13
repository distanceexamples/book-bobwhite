#' Function for integrating out the random effect
#' 
#' This function combines the poisson likelihoods for the observed counts given lambda and the 
#' normal density for the random effect coefficient
#' 
#' @param bj random effects coefficient
#' @param obs counts
#' @param lam lambda values before adding the random effect coefficient
#' @param std.ran random effect standard deviation

bl.function<-function(bl,obs,lam,std.ran){
  l<-length(bl)
  obs.prob<-array(NA,l)
  for (b in 1:l){
    lam2<-lam*exp(bl[b])
    obs.prob[b]<-prod(dpois(obs,lam2))*dnorm(bl[b],0,std.ran)
  }
  obs.prob  
}
