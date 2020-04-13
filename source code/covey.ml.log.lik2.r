#' Calculating the likelihood for model fitted to the covey data using ML methods
#' 
#' @details
#' The specific model this function is defined for is the preferred model from the RJMCMC algorithm
#' which has covariate Type in the detection model and covariates Type, JDctr and State in the count model
#' 
#' @param p parameters 
#' @param datax data object: covey.data
#' @param lim limits for integrating out the random effects

covey.ml.log.lik2<-function(p,datax,lim){
## Part 1 of covey.ml.log.lik2(): setting up the parameter values p for covariates 
# the detection model
scale.int<-p[1]     # the scale intercept 
shape<-exp(1.291197)# the fixed shape parameter is not estimated
sig.t<-c(0,p[2])    # coefficient for Type level "TREAT" (level "CONTROL is absorbed in the intercept)
# the count model 
int<-p[3]           # the intercept
typ<-c(0,p[4])      # coefficient for Type level "TREAT" (level "CONTROL" is absorbed in the intercept)
day<-p[5]           # coefficient for JDctr
st<-c(0,p[6:8])   # state coefficients for levels "MS","NC","TN" (level "MO" is absorbed in the intercept)
std.ran<-exp(p[9]) # the random effect standard deviation on the log-scale

# Part 2 of covey.ml.log.lik2(): the likelihood component pertaining to the detection model
# calculating the f(y) for each observed distances
le<-nrow(datax$dis.object)
fe<-numeric(le)
alltype<-sort(unique(datax$dis.object$Type))
dis.type<-match(datax$dis.object$Type,alltype)
# the sigma(z) for each detection
allscale<-exp(scale.int+sig.t[dis.type])
# calculating the f(y) for each observation (note that the truncation distance is stored in datax$w)
for (e in 1:le){
fe[e]<-f.haz.function.pt(datax$dis.object$distance[e],allscale[e],shape)/
       integrate(f.haz.function.pt,0,datax$w,allscale[e],shape)$value
}
# the sum of the log(f(y))
log.e<-sum(log(fe))

# Part 3 of covey.ml.log.lik2(): the likelihood component pertaining to the count model  
# setting up indices for the factor covariates and random effect coefficients
  State0<-match(datax$glmm.data$State,sort(unique(datax$glmm.data$State)))
  Type0<-match(datax$glmm.data$Type,sort(unique(datax$glmm.data$Type)))
  gr.id<-sort(unique(datax$glmm.data$gr.id))
  Ran0<-match(datax$glmm.data$gr.id,gr.id)
  J<-length(unique(datax$glmm.data$gr.id))
  # calculate the effective area for each visit to a point
  glmm.sig<-exp(scale.int+sig.t[Type0])#+sig.st[State0])
  n.ptvis<-nrow(covey.data$glmm.data) 
  l.efa<-array(NA,n.ptvis)
  for (j in 1:n.ptvis){
    l.efa[j]<-log(integrate(f.haz.function.pt,0,500,glmm.sig[j],shape)$value*pi*2)
  }

  ## calculate the log of the likelihood for each site while integrating out the random effect
  # the expected value for each count without the random effect coefficient
  lambda<-exp(int + (typ[Type0]) + (day*datax$glmm.data$JDctr) + st[State0]  + l.efa)
  lik<-array(NA,J)
  for (j in 1:J){
    j.rows<-which(Ran0==j)
    lik[j]<-integrate(bl.function,-lim,lim,obs=datax$glmm.data$detections[j.rows],lam=lambda[j.rows],std.ran=std.ran)$value
  }
  log.lik<-sum(log(lik))  
  log.lik.all<-log.lik+log.e
  # returns the negative log-likelihood as nlm minimizes
  -log.lik.all
}
