dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")
library(rlist)


input <- function(inputfile) {
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1]
  pfix = prefix()
  if (length(pfix) != 0) {
     prefix <- paste(pfix, "/", sep="")
  }
  clusterfile <- paste(pfix, toString(parameters["clusters", 2]), sep="")
  unthresholded <- paste(pfix, toString(parameters["unthresholded", 2]), sep="")
  thresholded <- paste(pfix, toString(parameters["thresholded", 2]), sep="")
  clusters <- read.csv(clusterfile, header=FALSE, check.names=FALSE);
  singlelist <<- clusters[,2]
  pos <- 0;
  listFull <<- list()
  for (i in 1:length(singlelist)) {
     if (is.na(clusters[,1][i])) {
        listFull <<- list.append(listFull, c());
        pos <- pos + 1;
     }
     else {
        listFull[[pos]] <<- append(listFull[[pos]], as.character(singlelist[[i]]));
     }
  }
  newdata <<- read.csv(thresholded, header=TRUE,  check.names=FALSE);
  originaldata <<- read.csv(unthresholded, header=TRUE,  check.names=FALSE);
  cn <<- colnames(newdata);
  #cn <<- cn[2:length(cn)];
  #newdata <<- newdata[,-1];
  newdata <<- apply(newdata, 1, as.numeric);
  newdata <<- t(newdata);
  colnames(newdata) <<- cn;
  cn <<- colnames(originaldata);
  #cn <<- cn[2:length(cn)];
  #originaldata <<- originaldata[,-1];
  originaldata <<- apply(originaldata, 1, as.numeric);
  originaldata <<- t(originaldata);
  colnames(originaldata) <<- cn;
}


run <- function() {
        mean_corr_full_inner_thresh <- c()
        mean_corr_full_inner_orig <- c()
        mean_corr_full_inner_thresh_neg0 <- c()
        mean_corr_full_inner_orig_neg0 <- c()
        mean_corr_full_cluster_matrix <<- matrix(nrow=length(listFull),ncol=length(listFull))
        std_corr_full_cluster_matrix <<- matrix(nrow=length(listFull),ncol=length(listFull))
        for(j in 1:length(listFull)){
                tmpBacteria <- listFull[[j]]
                tmpValues <- c()
                if(length(tmpBacteria)>1){
                        for(k in 1:(length(tmpBacteria)-1)){
                                for(m in (k+1):length(tmpBacteria)){
                                        cvalue <- newdata[tmpBacteria[k], tmpBacteria[m]]
                                        mean_corr_full_inner_thresh <- rbind(mean_corr_full_inner_thresh, c(tmpBacteria[k], tmpBacteria[m], abs(cvalue)));
                                        if(cvalue < 0)
                                                cvalue <- 0
                                        mean_corr_full_inner_thresh_neg0 <- rbind(mean_corr_full_inner_thresh_neg0, c(tmpBacteria[k], tmpBacteria[m], abs(cvalue)));
                                        cvalue <- originaldata[tmpBacteria[k], tmpBacteria[m]]
                                        tmpValues <- c(tmpValues, cvalue)
                                        mean_corr_full_inner_orig <- rbind(mean_corr_full_inner_orig, c(tmpBacteria[k], tmpBacteria[m], abs(cvalue)));
                                        if(cvalue < 0)
                                                cvalue <- 0
                                        mean_corr_full_inner_orig_neg0 <- rbind(mean_corr_full_inner_orig_neg0, c(tmpBacteria[k], tmpBacteria[m], abs(cvalue)));
                                }
                        }
                        mean_corr_full_cluster_matrix[j,j] <<- mean(tmpValues)
                        std_corr_full_cluster_matrix[j,j] <<- sd(tmpValues)
                }
        }

        mean_corr_full_outer_thresh <- c()
        mean_corr_full_outer_orig <- c()
        mean_corr_full_outer_thresh_neg0 <- c()
        mean_corr_full_outer_orig_neg0 <- c()
        for(j in 1:(length(listFull)-1)){
                tmpBacteriaJ <- listFull[[j]]
                for(jj in (j+1):length(listFull)){
                        tmpBacteriaJJ <- listFull[[jj]]
                        tmpValues <- c()
                        for(k in 1:length(tmpBacteriaJ)){
                                for(m in 1:length(tmpBacteriaJJ)){
                                        cvalue <- newdata[tmpBacteriaJ[k], tmpBacteriaJJ[m]]

                                        mean_corr_full_outer_thresh <- rbind(mean_corr_full_outer_thresh, c(tmpBacteriaJ[k], tmpBacteriaJJ[m], abs(cvalue)));
                                        if(cvalue < 0)
                                                cvalue <- 0
                                        mean_corr_full_outer_thresh_neg0 <- rbind(mean_corr_full_outer_thresh_neg0, c(tmpBacteriaJ[k], tmpBacteriaJJ[m], abs(cvalue)));
                                        cvalue <- originaldata[tmpBacteriaJ[k], tmpBacteriaJJ[m]]
                                        tmpValues <- c(tmpValues, cvalue)
                                        mean_corr_full_outer_orig <- rbind(mean_corr_full_outer_orig, c(tmpBacteriaJ[k], tmpBacteriaJJ[m], abs(cvalue)));
                                        if(cvalue < 0)
                                                cvalue <- 0
                                        mean_corr_full_outer_orig_neg0 <- rbind(mean_corr_full_outer_orig_neg0, c(tmpBacteriaJ[k], tmpBacteriaJJ[m], abs(cvalue)));
                                }
                        }
                        if(length(tmpValues)>0){
                                mean_corr_full_cluster_matrix[j,jj] <<- mean(tmpValues)
                                std_corr_full_cluster_matrix[j,jj] <<- sd(tmpValues)
                        }
                }
        }

        meansTable <<- c(
                mean(as.numeric(mean_corr_full_inner_thresh[,3])),
                sd(as.numeric(mean_corr_full_inner_thresh[,3])),
                mean(as.numeric(mean_corr_full_outer_thresh[,3])),
                sd(as.numeric(mean_corr_full_outer_thresh[,3])),

                mean(as.numeric(mean_corr_full_inner_orig[,3])),
                sd(as.numeric(mean_corr_full_inner_orig[,3])),
                mean(as.numeric(mean_corr_full_outer_orig[,3])),
                sd(as.numeric(mean_corr_full_outer_orig[,3])),

                mean(as.numeric(mean_corr_full_inner_thresh_neg0[,3])),
                sd(as.numeric(mean_corr_full_inner_thresh_neg0[,3])),
                mean(as.numeric(mean_corr_full_outer_thresh_neg0[,3])),
                sd(as.numeric(mean_corr_full_outer_thresh_neg0[,3])),

                mean(as.numeric(mean_corr_full_inner_orig_neg0[,3])),
                sd(as.numeric(mean_corr_full_inner_orig_neg0[,3])),
                mean(as.numeric(mean_corr_full_outer_orig_neg0[,3])),
                sd(as.numeric(mean_corr_full_outer_orig_neg0[,3]))
        )

}


output <- function(outputfile) {
        fileOut <- file.path(paste(outputfile, "meansFullMatrix.csv", sep="."));
        write.table(mean_corr_full_cluster_matrix, file=fileOut, sep=",", append=FALSE, col.names=NA, na="");
        fileOut <- file.path(paste(outputfile, "meStdevFullMatrix.csv", sep="."));
        write.table(std_corr_full_cluster_matrix, file=fileOut, sep=",", append=FALSE, col.names=NA, na="");
        fileOut <- file.path(paste(outputfile, "meansFull.csv", sep="."));
        write.table(meansTable, file=fileOut, sep=",", append=FALSE, col.names=NA, na="");

}

