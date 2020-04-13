#' Combined ikelihood for the full model fitted to the covey data
#' 
#' @details
#' This function is set up to fit the full model to the covey data where the detection model
#' contains the covariates Type and State and the count model contains the covariates Type, State and JDctr.
#' 
#' @param p Current model parameter values
#' @param raneff Current values for the random effect coefficients
#' @param datax Data object formatted with \code{create.data}: covey.data

covey.full.bayes.log.lik<-function(p, raneff, datax){
## Part 1 of covey.full.bayes.log.lik(): setting up the parameter values p for covariates 
# the detection model
scale.int<-p[1]     # the scale intercept 
shape<-exp(p[2])    # the shape parameter on the log scale
sig.t<-c(0,p[3])    # coefficient for Type level "Treat"" (level "CONTROL is absorbed in the intercept)
sig.st<-c(0,p[4:6]) # state coefficients for levels "MS","NC","TN" (level "MO" is absorbed in the intercept)
# the count model 
int<-p[7]           # the intercept
typ<-c(0,p[8])      # coefficient for Type level "Treat"" (level "CONTROL" is absorbed in the intercept)
day<-p[9]           # coefficient for JDctr
st<-c(0,p[10:12])   # state coefficients for levels "MS","NC","TN" (level "MO" is absorbed in the intercept)
std.ran<-exp(p[13]) # the random effect standard deviation on the log scale
# Part 2 of covey.full.bayes.log.lik(): the likelihood component pertaining to the detection model
# calculating the f(y) for each observed distances
le<-nrow(datax$dis.object)
fe<-numeric(le)
alltype<-sort(unique(datax$dis.object$Type))
dis.type<-match(datax$dis.object$Type,alltype)
allstate<-sort(unique(datax$dis.object$State))
dis.state<-match(datax$dis.object$State,allstate)
# the sigma(z) for each detection
allscale<-exp(scale.int+sig.t[dis.type]+sig.st[dis.state])
# calculating the f(y) for each observation (note that the truncation distance is stored in datax$w)
for (e in 1:le){
fe[e]<-f.haz.function.pt(datax$dis.object$distance[e],allscale[e],shape)/
       integrate(f.haz.function.pt,0,datax$w,allscale[e],shape)$value
}
# the sum of the log(f(y))
log.e<-sum(log(fe))

# Part 3 of covey.full.bayes.log.lik(): the likelihood component pertaining to the count model  
# setting up indices for the factor covariates and random effect coefficients
Type0<-match(datax$glmm.data$Type,sort(unique(datax$glmm.data$Type)))
State0<-match(datax$glmm.data$State,sort(unique(datax$glmm.data$State)))
gr.id<-sort(unique(datax$glmm.data$gr.id))
Ran0<-match(datax$glmm.data$gr.id,gr.id)
# calculate the effective area for each visit to a point
glmm.sig<-exp(scale.int+sig.t[Type0]+sig.st[State0])
n.ptvis<-nrow(covey.data$glmm.data) 
l.efa<-array(NA,n.ptvis)
for (j in 1:n.ptvis){
  l.efa[j]<-log(integrate(f.haz.function.pt,0,datax$w,glmm.sig[j],shape)$value*pi*2)
}
# calculate the log of the Poisson likelihood for each count
lambda<-exp(int + typ[Type0] + (day*datax$glmm.data$JDctr) + st[State0] + raneff[Ran0] + l.efa)
dpois.y<-dpois(datax$glmm.data$detections,lambda)
logdpois.y<-sum(log(dpois.y))
# calculate the log of the normal density for each random effect coefficient
log.dnorm.raneff<-sum(log(dnorm(raneff,0,std.ran)))
# adding up the likelihood components
log.lik<-logdpois.y + log.dnorm.raneff + log.e
# the function return
return(log.lik)
}

