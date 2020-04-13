#' Updating the random effect coefficients
#' 
#' @details
#' This function is based on covey.full.bayes.log.lik and updates all random effect coefficients. 
#' However, it does not calculate the likelihood for the detection model. 
#' It only calculates the offset once and incorporates it into the lambda model. 
#' In contrast to covey.full.bayes.log.lik which returns the log-likelihood, this function returns
#' the updated random effect coefficients. 
#' 
#' @param p Current model parameter values
#' @param raneff Current values for the random effect coefficients
#' @param datax Data object formatted with \code{create.data}: covey.data

covey.full.bayes.log.lik.raneff.update<-function(p,raneff,datax){
  # setting up the covariates
  scale.int<-p[1]
  shape<-exp(p[2])
  sig.t<-c(0,p[3]) # coefficient for type = treat
  sig.st<-c(0,p[4:6])
  # 4 x 2 matrices: row = state, columns 1:2 type = control, treat
  sig.t.mat<-matrix(c(rep(0,1),rep(p[3],1)),4,2,byrow=T)
  sig.st.mat<-matrix(c(0,p[4:6]),4,2,byrow=F)
  int<-p[7]
  typ<-c(0,p[8])
  day<-p[9]
  st<-c(0,p[10:12])
  std.ran<-exp(p[13])
  
#   # calculating the f(y) for each observed distances
#   le<-nrow(datax$dis.object)
#   fe<-numeric(le)
#   alltype<-sort(unique(datax$dis.object$Type))
#   dis.type<-match(datax$dis.object$Type,alltype)
#   allstate<-sort(unique(datax$dis.object$State))
#   dis.state<-match(datax$dis.object$State,allstate)
#   allscale<-exp(scale.int+sig.t[dis.type]+sig.st[dis.state])
#   # calculating the pdf f(y) for each observation
#   for (e in 1:le){
#     fe[e]<-f.haz.function.pt(datax$dis.object$distance[e],allscale[e],shape)/integrate(f.haz.function,0,500,allscale[e],shape)$value
#   }
#   # the sum of the log(f(y))
#   log.e<-sum(log(fe))

n.ptvis<-nrow(covey.data$glmm.data) 
l.efa<-array(NA,n.ptvis)
State0<-match(datax$glmm.data$State,sort(unique(datax$glmm.data$State)))
Type0<-match(datax$glmm.data$Type,sort(unique(datax$glmm.data$Type)))
gr.id<-sort(unique(datax$glmm.data$gr.id))
Ran0<-match(datax$glmm.data$gr.id,gr.id)
glmm.sig<-exp(scale.int+sig.t[Type0]+sig.st[State0])
for (j in 1:n.ptvis){
  l.efa[j]<-log(integrate(f.haz.function.pt,0,datax$w,glmm.sig[j],shape)$value*pi*2)
}

lambda<-exp(int + (typ[Type0]) + (day*datax$glmm.data$JDctr) + st[State0] + raneff[Ran0] + l.efa)
#lambda<-as.numeric(unlist(lambda))
for (r in 1:length(gr.id)){
 u<-rnorm(1,0,0.2)
 yr<-which(datax$glmm.data$gr.id==gr.id[r])
 dpois.y.new<-dpois(datax$glmm.data$detections[yr],lambda[yr]*exp(u))
 num<-sum(log(dpois.y.new))+log(dnorm(raneff[r]+u,0,std.ran))
 dpois.y.cur<-dpois(datax$glmm.data$detections[yr],lambda[yr])
 den<-sum(log(dpois.y.cur))+log(dnorm(raneff[r],0,std.ran))
    A<-min(1,exp(num-den))
    V<-runif(1)
    if(V<=A){raneff[r]<-raneff[r]+u}
}
raneff
}
  
  
  