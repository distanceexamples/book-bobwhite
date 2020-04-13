#' Format data
#'
#' This function creates data frames with the format and columns needed for integrated likelihood modeling with random effect count models
#' The data frames are \code{dis.object} which is at the resolution of the detections and \code{glmm.data} which is at the resolution of the counts at the sampler (or counts in the distance bins).

#' @param data data frame from user; 
#' @param w truncation distance in units specified in \code{dis.unit}
#' @param sampler "lines" (default) or "points" 
#' @param dis.unit units of distances (perpendicular for lines and radial for points)
#' @param sampler.units required field for line transects (default = \code{dis.unit})
#' @param exact TRUE for exact distance measurements, FALSE for interval distance data
#' @param cutpoints required if exact = FALSE: list of cutpoints in same units as \code{dis.unit}, needs to be of length I+1 for I cutpoints with cutpoints[1] = 0 and cutpoints[I+1] = w
#' @param dis.cov optional list of covariates that may be included in MCDS model (expects a list of column names from \code{data} where covariates are found)
#' @param count.cov required list of covariates that may be included in count model (expects the column names from data where covariates are found, no NA's allowed, only one value for each visit to a sampler allowed)
#' @param dis.factor.cov list of covariates to be considered as factors
#' @param count.factor.cov list of covariates to be considered as factors
#' @param conditional TRUE creates 1 column of for detections and 1 column for individuals in glmm.data. If False, creates 1 column for detections and individuals for each distance bin (only applicable for interval distance data)

#' @details 
#' 
#' For \code{data} the following columns are expected:
#' 
#'  \code{$gr.id} identifiers for groups (grouping factor for random effect); 
#'  
#'  \code{$smp.id} identifiers for points/lines; 
#'  
#'  \code{$vis.id} identifiers for individual visits to the same sampler (in case of repeat visits); different visits need to have different lables for each unique smp.id;
#'  
#'  \code{$det.id} identifiers for detections; 
#'  
#'  \code{$line.length} for line transects only; 
#'  
#'  \code{$distance} exact or midpoints for intervals for interval data; 
#'  
#'  \code{$size} cluster size; 
#'  
#'  covariates for detection or count model (no NA's allowed);
#'  
#' Also for \code{data} the following rows are expected: 
#' 
#' 1 row for each detection; 
#' 
#' if no detection at a visit to a sampler: 1 row with NA's for \code{$det.id}, \code{$distance} and \code{$size}

#' @examples
#' covey.data<-create.data(covey,500,"points","m","m",F,dis.cov=c("Type","State"),count.cov=c("Type","State","JDctr"), dis.factor.cov = c("Type","State"), count.factor.cov = c("Type","State"), conditional=T)

create.data <-
function(data, w, sampler = "lines", dis.unit, sampler.unit = dis.unit, 
  binned, cutpoints = NULL, dis.cov = NULL, count.cov, dis.factor.cov = NULL,
  count.factor.cov = NULL, conditional = T) {
        
        data$smp.vis<-paste(data$gr.id,data$smp.id,data$vis.id,sep="_") # creates a unique identifier for each visit to each sampler of each group
        
        data.object<-NULL                                               # object that will be returned by this function
        data.object$w <- w                                              # the truncation distance
        data.object$sampler <- sampler                                  # lines or points
        data.object$dis.unit <- dis.unit                                # units for distances (cm, m, km)
        data.object$sampler.unit <- sampler.unit                        # units for $line.length (cm, m, km) 
        data.object$binned <- binned
        data.object$cutpoints <- cutpoints
        
        data.object$dis.cov <- dis.cov
        dis.cov.fac <- array(NA,length(dis.cov))
        for (u in 1:length(dis.cov)){
          dis.cov.fac[u] <- ifelse(is.na(match(dis.cov[u],dis.factor.cov))==F,TRUE,FALSE)
        }
        data.object$dis.cov.fac<-dis.cov.fac
        
        data.object$count.cov<-count.cov
        count.cov.fac<-array(NA,length(count.cov))
        for (u in 1:length(count.cov)){
          count.cov.fac[u]<-ifelse(is.na(match(count.cov[u],count.factor.cov))==F,TRUE,FALSE)
        }
        data.object$count.cov.fac<-count.cov.fac
        
        data.object$conditional<-conditional
        
        # $dis.object = a data.frame containing information required for detection function model L_{y|z}
        if(sampler == "lines")  {dis.object<-data[which(is.na(data$det.id)==F),c("smp.vis","det.id","smp.id","gr.id","vis.id","line.length","distance","size",dis.cov)]}
        if(sampler == "points") {dis.object<-data[which(is.na(data$det.id)==F),c("smp.vis","det.id","smp.id","gr.id","vis.id","distance","size",dis.cov)]  }
        # truncating at w
          if(binned==T){
            bins<-which.bin(data$distance,cutpoints)
            data$bin<-bins[[1]]
            data[,c("distbegin","distend")]<-bins[[2]]
        }
        dis.object<-data[which(data$distance<=w),]
 
       # adding to the data.object
        data.object$dis.object<-data.frame(dis.object)
  
       # $glmm.data a data frame containing the information for the random effect count model L_n
       if(conditional==T){
        data.object$glmm.data<-create.glmm.data(data, w, sampler, count.cov, count.cov.fac) }
        else{
        data.object$glmm.data<-create.glmm.data.unconditional(data, w, sampler, count.cov, count.cov.fac)           
        }
        data.object
    }
          
          

