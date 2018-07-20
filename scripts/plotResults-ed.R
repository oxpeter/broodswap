
# Function to create plots during optimization process.

plotResults <- function(result, data, title){
  print("Clust Size Max Min Diff")
  plot.fun <- function(x)    lines(x=x$time, y=x$y)
  temp <- map(result$z)
  for(j in unique(temp)){
    comp.ids <- which(temp == j)
    comp.data <- na.omit(data[which(data$indiv %in% comp.ids),])

    # Check that there is more than one curve in the cluster and then plot mean plus curves
    if(length(comp.ids) > 1) {
      write.table(as.data.frame(unique(comp.data$gene)), file=paste(title, ".cluster_", j, ".list", sep=""), quote=FALSE)
      pdf(paste(title, ".cluster_",j,".pdf", sep=""))
      init.data <- subset(comp.data, indiv==comp.ids[1])

      plot(x=init.data$time, y=init.data$y, type="n", ylim=range(comp.data$y),
           ylab="Log2 (Fold Expression Level)",
           xlab="Time",
           main=paste("Cluster ", j, " (", length(comp.ids), " genes)",  sep=""))

      # plot spline curve
      by(comp.data, comp.data$indiv, plot.fun)
      #est.fit <- getF(result$fit[[j]])$f.time[[1]]
      #lines(est.fit[,1], est.fit[,2], col="orange", lwd=3)

      # calculate mean expression and plot points
      mean.expression <- tapply(comp.data$y,comp.data$time, mean)
      points(x=init.data$time,
             y=mean.expression,
             col="red",
             pch=19
             )

      dev.off()

      max_y <- max(mean.expression)
      min_y <- min(mean.expression)
      print(sprintf("%-3d %-4d %-5.2f %-5.2f %-5.2f",
                    j, length(comp.ids),max_y,min_y,max_y - min_y))
    }

    # If only one curve in a cluster, just plot the curve - no mean.
    else{
      write.table(as.data.frame(unique(comp.data$gene)), file=paste(title, ".cluster_", j, ".list", sep=""), quote=FALSE)
      pdf(paste(title, ".cluster_",j,".pdf", sep=""))
      comp.data <- na.omit(data[which(data$indiv %in% comp.ids),])
      plot(x=comp.data$time, y=comp.data$y,
           type="l", ylab="Expression level", xlab="Time", main=paste("Cluster ", j, " (", length(comp.ids), " genes)",  sep=""))
      dev.off()
    }
  # functions to pause graph until user is ready to display the next cluster:
  #Sys.sleep(5)
  #cat ("Press [enter] to continue")
  #line <- readline()
  }
}

pltResultsTest <- function(result, data, title){
    print("Test works. Good luck on the rest.")
    }