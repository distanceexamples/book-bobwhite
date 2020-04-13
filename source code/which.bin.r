#' Determining the distance bin 
#'
#' This function takes a single distance and determines which distance bin it belongs to
#'
#' @param x vector of perpendicular distances from the line or radial distances from the point
#' @param cutpoints a vector of cutpoint of the intervals
#'
#' @details
#' if a value in x matches a cutpoint in cutpoints the value of x will be attributed to the bin closer to the line/point unless the value in x is 0
#'
#' @value
#' list of three objects: 
#' 1. a matrix of dimensions I x 2, where I equals the number of bins, containing the cutpoints for each bin. 
#' 2. a vector of length(x) containing the number of which bin the respective distance in vector \code{x} belongs to 
#' 3. a matrix of dimensions length(x) x 2 containing the cutpoints of the bin that the distances in \code{x} belong to (the values needed for \code{distbegin} and \code{distend}

which.bin<-function(x,cutpoints){
bins<-matrix(0,length(cutpoints)-1,2)
for (i in 1:(length(cutpoints)-1)){
bins[i,1:2]<-cutpoints[i:(i+1)]
}
y<-array(NA,length(x))
yz<-matrix(NA,length(x),2)
for (i in which(is.na(x)==F)){
      if (x[i]==0) {y[i]<-1}
      else {if(x[i]<=max(cutpoints)){
          y[i]<-which(bins[,1]<x[i] & bins[,2]>=x[i])
          yz[i,]<-bins[y[i],]}}
      }
return(list(y,yz,bins))
}

