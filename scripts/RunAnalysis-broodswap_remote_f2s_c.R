source('plotResults-ed.R')
source('EMAlgPsplineLargeData.R')
# Packages required.
require(mclust)
require(amer)
require(DAAG)
require(mvtnorm)


# Set up a dataframe using a broodswap dataset.

# times are log transformed from c(0,6,12,24,48,96)   (xt = ln(x+1))
gene.time <- c(0,2.56494936,3.21887582,3.8918203,4.57471098)
indiv.time <- c(0,2.56494936,3.21887582,3.8918203,4.57471098)

full.data <- read.table("../create_data_tables/forg_norm_vst_idx.tbl", header=T)
full.data <- data.frame(y=c(t(full.data[,3:7])),
                        time=rep(indiv.time, nrow(full.data)),
                        indiv=factor(rep(1:nrow(full.data),each=5)),
                        gene=rep(full.data[,2],each=5)
                        )

# Set default values for clustering algorithm.
maxiters <- 15  # number of iterations for RCEM step
thresh <- 0.5   # Threshold (c) to use in RCM step

# Number of clusters to fit
n.comps <- c(42,52,44,54,46,56,48,58,40,50)

est <- "REML"   # Use REML to fit the mixed model.
subset.idx <- 1:length(unique(full.data$indiv))

# Function to choose the number and location of the knots for use with amer (instead of default values).
default.knots <- function(x,num.knots)
{
   if (missing(num.knots))
      num.knots <- max(5,min(floor(length(unique(x))/4),35))
   return(quantile(unique(x),seq(0,1,length=(num.knots+2))[-c(1,(num.knots+2))], na.rm=T))
}

# Determine the number and location of knots.
my.knots <- default.knots(full.data$time)

# Function to create an initial random start for the clustering algorithm. Randomly allocate
# each gene to ncomps clusters.
createRandomStart <- function(ncomps, subset.idx, gene.time){
  rstart<-sample(1:ncomps,length(subset.idx),replace=T)
  r<-1
  while(any(table(rstart)==1)==T && r<100){
    r<-r+1
    rstart<-sample(1:ncomps,length(subset.idx),replace=T)
  } # end while
  temp <- rep(rstart, each=length(gene.time))
  z <- unmap(temp)
  return(z)
}

# Store the results.
result.sim <- list()
temp.sim <- list() # Store the timings. (May not be required.)
time.sim <- c()
result.time <- c()


# For each cluster size given in n.comps, fit the model and generate the clusters.
for(ncomps in n.comps){
    # Take 5 random starts and store the results from the EM algorithm. Attempts to avoid local minimum.
    for(n.starts in 1:5){

      print(paste("n.components = ", ncomps, " n.starts = ", n.starts, sep=""))
      test_random <- createRandomStart(ncomps, subset.idx, gene.time)

      # Generate a random start.
      z <- Matrix(createRandomStart(ncomps, subset.idx, gene.time),sparse=T)

      # Try to cluster and store fit and CPU time.
      time.sim[n.starts] <- system.time(temp.sim[[n.starts]] <- try(EMAlg.pspline.large(data=full.data, EMmaxit = maxiters, ncomps = ncomps, z=z, random="1|indiv", est=est, n.time=length(gene.time), knots=my.knots, thresh=thresh)))[[3]]


      restarts <- 0
      # If there's an error during fitting try 2 more restarts before recording
      # results as NA
      while(is.character(temp.sim[[n.starts]][1])&& substr(temp.sim[[n.starts]][1],1,5) == "Error" && restarts <= 2){
        z <- Matrix(createRandomStart(ncomps, subset.idx, gene.time),sparse=T)

        time.sim[n.starts] <- system.time(temp.sim[[n.starts]] <- try(EMAlg.pspline.large(data=full.data, EMmaxit = maxiters, ncomps = ncomps, z=z, random="1|indiv", est=est, n.time=length(gene.time),  knots=my.knots, thresh=thresh)))[[3]]
        restarts <- restarts + 1
      }

      if(is.character(temp.sim[[n.starts]][1])&& substr(temp.sim[[n.starts]][1],1,5) == "Error" && restarts > 2){
          time.sim[n.starts] <- NA
          temp.sim[[n.starts]] <- NA
      }

      # Pick the best solution and corresponding time from 5 random starts.
      temp.sim <- temp.sim[which(!is.na(temp.sim))]
      like<-lapply(temp.sim, function(x) x$L)
      result.sim <- temp.sim[[which(unlist(like) == max(unlist(like),na.rm=T))]]#[c(2,3,5)]
      result.time <- sum(time.sim, na.rm=T)
    }
    result <- list(result.sim, result.time)
    save(result, file=paste("../som_forg2stat/result", ncomps, "comps.f2s_logtime.RData", sep="")) # Save the result for each number of clusters	  # specified in n.comps.
}

## Plot the BIC curve to choose the number of clusters.
#BIC <- c()
#for(i in 1:length(n.comps)){
#    load(paste("../som_forg2stat/result", n.comps[i], "comps.f2s_logtime.RData", sep=""))
#    BIC[i] <- result[[1]]$BIC
#    print(paste(n.comps[i], "-->", BIC[i], sep=" "))
#}
#plot(n.comps, BIC, type="b", xlab="G", ylab="BIC")
## Plot the clusters for the best result (as chosen by BIC)
## For example, choose the solution for 6 clusters.
#
#best_G <- 110
#best_file <- paste("../som/result", best_G, "comps.f2s_logtime.RData", sep="")
#load(best_file)
#plotResults(result[[1]], data=full.data, "results_110_clust.S2F")