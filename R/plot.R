#' @title Plot the clustering results of subgroups. 
#' @description plot the observed response rate versus subgroup ID with clusters coded by the color of dots. 
#' @param res BCHM calculation results.
#' @param col Color vector
#' @param pch pch vector
#' @param cex size of points

#' @examples
#' nDat = c(25, 25, 25, 25) # total number of patients
#' xDat = c(2, 3, 8, 6)  # number of responses
#' alpha <- 1e-20 
#' d0 <- 0.0 
#' alpha1 = 50   
#' beta1 = 10  
#' tau2 <- 0.1  
#' phi1 <- 0.1  
#' deltaT <- 0.2  
#' thetaT <- 0.60   
#' 
#' res <- BCHM(nDat = nDat,
#'             xDat = xDat,
#'             alpha = alpha,
#'             d0 = d0,             
#'             alpha1 = alpha1, 
#'             beta1 = beta1,
#'             tau2 = tau2,
#'             phi1 = phi1, 
#'             deltaT = deltaT,
#'             thetaT = thetaT,
#'             burnIn = 100,
#'             MCIter = 200,
#'             MCNum = 1000,
#'             seed = 1000
#' )
#' print(res$SMatrix)
#' print(res$Result)
#' col <- res$Result[,4]
#' 
#' BCHMplot_cluster(res, col, pch=16)
#'
#' @return None 
#' @seealso \code{\link{BCHM} Perform the analysis based on the BCHM design.}
#' @seealso \code{\link{BCHMplot_post_value} Plot the posterior response of subgroups. }
#' @seealso \code{\link{BCHMplot_post_dist} Plot the posterior distributions of subgroups. }

#' @export
BCHMplot_cluster <- function(res, col, pch, cex=2){
  UseMethod("BCHMplot_cluster")
}

#' @importFrom crayon red
#' @export
BCHMplot_cluster.default <- function(res, col, pch, cex=2) {
  stop(red(
    "Don't know how to make a plot with an object of type",
    paste(class(res), collapse = ", "), "."))
}



#' @importFrom stats median var sd runif density rnorm
#' @importFrom graphics plot
#' @importFrom graphics axis legend lines points
#' 
#' @export
BCHMplot_cluster.BCHM_result <- function(res, col, pch, cex=2)
{
  s <- res$Result
  
  plot(
    1:dim(s)[1],
    s[, 3],
    xlab="Subgroup ID", ylab="Observed Response Rates",
    main="Subgroup Clusters",
    xlim = c(0.5, dim(s)[1] + 0.5),
    pch = pch,  #s[, 4] + 15,
    col = col, #s[, 4],
    cex = cex #2 #weight / 10
  )
  xtick<- 1:dim(s)[1]
  axis(side=1, at=xtick)
}



#' @title Plot the posterior response of subgroups. 
#' @description plot the posterior response rate with its highest probability density (HPD) interval by subgroup ID
#' @param res BCHM calculation results.
#' @param col Color vector
#' @param pch pch vector pch[1] Posterior mean pch[2] Observed mean
#' @param cex size of points
#' @param HPD Highest Posterior Density level for drawing (NA: No HPD drawing )
#' @param ObsMean Draw the observed mean of subgroups if this parameter is TRUE
#' @examples
#' nDat = c(25, 25, 25, 25) # total number of patients
#' xDat = c(2, 3, 8, 6)  # number of responses
#' alpha <- 1e-20 
#' d0 <- 0.0 
#' alpha1 = 50   
#' beta1 = 10  
#' tau2 <- 0.1  
#' phi1 <- 0.1  
#' deltaT <- 0.2  
#' thetaT <- 0.60   
#' 
#' res <- BCHM(nDat = nDat,
#'             xDat = xDat,
#'             alpha = alpha,
#'             d0 = d0,             
#'             alpha1 = alpha1, 
#'             beta1 = beta1,
#'             tau2 = tau2,
#'             phi1 = phi1, 
#'             deltaT = deltaT,
#'             thetaT = thetaT,
#'             burnIn = 100,
#'             MCIter = 200,
#'             MCNum = 1000,
#'             seed = 1000
#' )
#' print(res$SMatrix)
#' print(res$Result)
#' col <- res$Result[,4]
#' 
#' BCHMplot_post_value(res, col, HPD = 0.8)
#' 
#' @return None
#' @seealso \code{\link{BCHM} Perform the analysis based on the BCHM design.}
#' @seealso \code{\link{BCHMplot_cluster} Plot the clustering results of subgroups. }
#' @seealso \code{\link{BCHMplot_post_dist} Plot the posterior distributions of subgroups. }
#' @export
BCHMplot_post_value <- function(res, col, pch=c(19, 4), cex=2, HPD = 0.95, ObsMean = FALSE) {
  UseMethod("BCHMplot_post_value")
}

#' @importFrom crayon red
#' @export
BCHMplot_post_value.default <- function(res, col, pch=c(19, 4), cex=2, HPD = 0.95, ObsMean = FALSE){
  stop(red(
    "Don't know how to make a plot with an object of type",
    paste(class(res), collapse = ", "), "."))
}


#' @importFrom stats median var sd runif density rnorm
#' @importFrom graphics plot
#' @importFrom graphics axis legend lines points
#' 
#' @export
BCHMplot_post_value.BCHM_result <- function(res, col, pch=c(19, 4), cex=2, HPD = 0.95, ObsMean = FALSE)
{
  if (length(cex) < 2)
  {
    cex <- rep(cex, 2)
  }
  s <- res$Result
  if (is.na(HPD))
  {
    plot(
      1:dim(s)[1],
      s[, 5],
      xlab = "Subgroup ID",
      ylab = "Posterior Response Rates",
      main = "Posterior Probability",
      xlim = c(0.5, dim(s)[1] + 0.5),
      ylim = c(0, min(1, max(s[, 5]) + 0.1)),
      pch = pch[1], #19 s[, 4] + 15,
      col = col,
      #s[, 4],
      cex = cex[1] #weight / 10
    )
    xtick<- 1:dim(s)[1]
    axis(side=1, at=xtick)
  } else {
    
    plot(
      1:dim(s)[1],
      s[, 5],
      xlab = "Subgroup ID",
      ylab = "Posterior Response Rates",
      main = paste0("Posterior Probability HPD = ", HPD),
      xlim = c(0.5, dim(s)[1] + 0.5),
      ylim = c(0, min(1, max(s[, 5]) + 0.3)),
      pch = pch[1], #s[, 4] + 15,
      col = col, #s[, 4],
      cex = cex[1] #weight / 10
    )
    xtick<- 1:dim(s)[1]
    axis(side=1, at=xtick)
    samp <- res$Samples
    for (i in 1:dim(samp)[2])
    {
      hpdLevel <- boa.hpd(samp[,i], 1 - HPD)
      lines(c(i, i), c(hpdLevel$lower, hpdLevel$upper), col = col[i], lwd = 3)
      wd <- 0.1
      lines(c(i - wd, i + wd), c(hpdLevel$lower, hpdLevel$lower), col = col[i], lwd = 3)
      lines(c(i - wd, i + wd), c(hpdLevel$upper, hpdLevel$upper), col = col[i], lwd = 3)      
      #browser()
    }
  }
  if(ObsMean)
  {
    points(
      1:dim(s)[1],
      s[, 3],
      pch = pch[2], #4, #19, #s[, 4] + 15,
      col = "magenta",#col, #s[, 4],
      cex = cex[2] #weight / 10
    )
  }
  
}



#' @title Plot the posterior distributions of subgroups. 
#' @description plot the posterior distribution by subgroup ID
#' @param res BCHM calculation results.
#' @param col Color vector
#' @param lty line types
#' @param lwd line width
#' @param xlim X-axis range
#' @param ylim Y-axis range
#' @examples
#' nDat = c(25, 25, 25, 25) # total number of patients
#' xDat = c(2, 3, 8, 6)  # number of responses
#' alpha <- 1e-20 
#' d0 <- 0.0 
#' alpha1 = 50   
#' beta1 = 10  
#' tau2 <- 0.1  
#' phi1 <- 0.1  
#' deltaT <- 0.2  
#' thetaT <- 0.60   
#' 
#' res <- BCHM(nDat = nDat,
#'             xDat = xDat,
#'             alpha = alpha,
#'             d0 = d0,             
#'             alpha1 = alpha1, 
#'             beta1 = beta1,
#'             tau2 = tau2,
#'             phi1 = phi1, 
#'             deltaT = deltaT,
#'             thetaT = thetaT,
#'             burnIn = 100,
#'             MCIter = 200,
#'             MCNum = 1000,
#'             seed = 1000
#' )
#' print(res$SMatrix)
#' print(res$Result)
#' col <- res$Result[,4]
#' 
#' BCHMplot_post_dist(res, col, lty=1:length(nDat), lwd =3, xlim=c(0, 0.8))
#' 
#' @return None
#' @seealso \code{\link{BCHM} Perform the analysis based on the BCHM design.}
#' @seealso \code{\link{BCHMplot_cluster} Plot the clustering results of subgroups. }
#' @seealso \code{\link{BCHMplot_post_value} Plot the posterior response of subgroups. }
#' @export
BCHMplot_post_dist <- function(res, col, lty, lwd = 2, xlim=c(0,1), ylim = NA){
  UseMethod("BCHMplot_post_dist")
}



#' @importFrom crayon red
#' @export
BCHMplot_post_dist.default <- function(res, col, lty, lwd = 2, xlim=c(0,1), ylim = NA){
  stop(red(
    "Don't know how to make a plot with an object of type",
    paste(class(res), collapse = ", "), "."))
}


#' @importFrom stats median var sd runif density rnorm
#' @importFrom graphics plot
#' @importFrom graphics axis legend lines points
#' @export
BCHMplot_post_dist.BCHM_result <- function(res, col, lty, lwd = 2, xlim=c(0,1), ylim = NA)
{
  r <- res$Result
  #xLim <- max(r[, 2] / r[, 1]) + 0.3
  s <- res$Samples
  
  maxY <- 0
  for(i in 1:dim(s)[2])
  {
    sampledP <- s[,i]
    maxY <- max(maxY, max(density(sampledP)$y))
  }
  if (sum(is.na(ylim)) > 0)
  {
    ylim <- c(0, maxY * 1.1)
  }
  legendStr <- c()
   plot(c(0,0), ylim = ylim, xlim = xlim, col = "white", xlab="Response Rates", xaxt='n', ylab="Density", main="Posterior Distribution")
   xtick<- (0:5) / 5
   axis(side=1, at=xtick)   
   for (i in 1:dim(s)[2])
   {
     lines(density(s[,i]),
           col = col[i], lty=lty[i], lwd = lwd)
     legendStr <- c(legendStr, paste("Subg.", i))
   }

   legend("topright", legendStr, lty=lty, lwd = lwd, col = col)
}

