#' Matching strings of numbers against a matrix
#' 
#' This function matches a string of numbers against a matrix
#' and returns the row number of the matrix that matches the string exactly
#' 
#' @param xmod a string of numbers
#' @param model.matrix A matrix of numbers

match.function<-function(xmod,model.matrix){
nmodel<-nrow(model.matrix)
npar<-ncol(model.matrix)
is.match.or.not<-array(NA,nmodel)
for (i in 1:nmodel){
is.match.or.not[i]<-sum(xmod==model.matrix[i,])}
result<-which(is.match.or.not==npar)
return(result)
}
