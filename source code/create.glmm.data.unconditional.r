#' Format distance data into count data
#'
#' This function converts distance sampling data containing records of detections, into count data, where each visit to a sampler are the records.
#' It is an internal function that is called by the \code{creat.data} function.
#'
#' @param data data frame from user;
#' expects the following columns:
#'  \code{$gr.id} (identifiers for groups (grouping factor for random effect));
#'  \code{$smp.id} (identifiers for points/lines);
#'  \code{$vis.id} (identifiers for individual visits to the same sampler (in case of repeat visits));
#'  \code{$det.id} (identifiers for detections);
#'  \code{$line.length} (for line transects only);
#'  \code{$distance} (exact or midpoints for intervals for interval data);
#'  \code{$size} (cluster size);
#'  covariates for count model (no NA's allowed);
#' expects the following rows:
#'  1 row for each detection;
#'  if no detection at a visit to a sampler: 1 row with NA's for \code{$det.id}, \code{$distance} and \code{$size}
#' @param w truncation distance in units specified in \code{dis.unit}
#' @param sampler "lines" (default) or "points"
#' @param count.cov required list of covariates that may be included in count model (expects the column names from \code{data}, no NA's allowed, only one value for each visit to a sampler allowed)
#' @param count.cov.fac list of covariates in \code{count.cov} that are factors

create.glmm.data.unconditional<-function(data, w, sampler, count.cov, count.cov.fac) {
    data$smp.vis<-paste(data$gr.id,data$smp.id,data$vis.id,sep="_") # creates a unique identifier for each visit to each sampler of each group
    glmm.raw<-data[,c("smp.vis","det.id","smp.id","gr.id","vis.id","distance","size",count.cov)]
    if(sampler=="lines"){glmm.raw$line.length<-data[,"line.length"]}
    # truncating the data at w while making sure that effort lines remain in place
    truncated<-which(glmm.raw$distance>w)
    if(length(truncated)>0){
        glmm.raw[truncated,c("det.id","distance","size")]<-NA
    }
    # list of unique smp.vis (each unique smp.vis will have one row in glmm.data)
    smp.vis<-unique(data$smp.vis)

    if(sampler=="lines"){glmm.data<-matrix(NA,length(smp.vis),7+length(count.cov))
    glmm.data<-data.frame(glmm.data)
    glmm.data$smp.vis<-smp.vis
    for (s in 1:length(smp.vis)){
        rows<-which(glmm.raw$smp.vis==smp.vis[s])
        glmm.data$smp.id[s]<-paste(glmm.raw$smp.id[rows[1]])
        glmm.data$gr.id[s]<-paste(glmm.raw$gr.id[rows[1]])
        glmm.data$vis.id[s]<-paste(glmm.raw$vis.id[rows[1]])
        glmm.data$detections[s]<-length(which(is.na(glmm.raw$det.id[rows])==F))
        glmm.data$individuals[s]<-sum(glmm.raw$size[rows],na.rm=T)
        glmm.data$line.length[s]<-paste(glmm.raw$line.length[rows[1]])
        for (t in 1:length(count.cov)){
          if(count.cov.fac[t]==T){
            glmm.data[s,count.cov[t]]<-paste(glmm.raw[rows[1],count.cov[t]])
          }
          else{
            glmm.data[s,count.cov[t]]<-glmm.raw[rows[1],count.cov[t]]
          }
        }
    }
    }
    
    if(sampler=="points"){glmm.data<-matrix(NA,length(smp.vis),6+length(count.cov))
    glmm.data<-data.frame(glmm.data)
    colnames(glmm.data)<-c("smp.vis","smp.id","gr.id","vis.id","detections","individuals",count.cov)
    glmm.data$smp.vis<-smp.vis
    for (s in 1:length(smp.vis)){
        rows<-which(glmm.raw$smp.vis==smp.vis[s])
        glmm.data$smp.id[s]<-paste(glmm.raw$smp.id[rows[1]])
        glmm.data$gr.id[s]<-paste(glmm.raw$gr.id[rows[1]])
        glmm.data$vis.id[s]<-paste(glmm.raw$vis.id[rows[1]])
        glmm.data$detections[s]<-length(which(is.na(glmm.raw$det.id[rows])==F))
        glmm.data$individuals[s]<-sum(glmm.raw$size[rows],na.rm=T)
        for (t in 1:length(count.cov)){
          if(count.cov.fac[t]==T){
            glmm.data[s,count.cov[t]]<-paste(glmm.raw[rows[1],count.cov[t]])
          }
          else{
            glmm.data[s,count.cov[t]]<-glmm.raw[rows[1],count.cov[t]]
          }
        }
    }
    }
  for (i in 1:(length(cutpoints)-1)){
        detname<-paste("det",i,sep="_")
        nname<-paste("n",i,sep="_")
  for (s in 1:length(smp.vis)){
        rows<-which(glmm.raw$smp.vis==smp.vis[s] & glmm.raw$bin==i)
        glmm.data[s,detname]<-length(which(is.na(glmm.raw$det.id[rows])==F))
        glmm.data[s,nname]<-sum(glmm.raw$size[rows],na.rm=T)
  }
  }
    return(glmm.data)
}


