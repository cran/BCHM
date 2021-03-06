% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{BCHMplot_post_value}
\alias{BCHMplot_post_value}
\title{Plot the posterior response of subgroups.}
\usage{
BCHMplot_post_value(
  res,
  col = res$Result[, 4],
  pch = c(19, 4),
  cex = 2,
  HPD = 0.95,
  xlim = c(0, dim(res$Result)[1] + 2),
  ylim = c(0, 1),
  ObsMean = FALSE,
  ...
)
}
\arguments{
\item{res}{BCHM calculation results.}

\item{col}{Color vector}

\item{pch}{pch vector pch[1] Posterior mean pch[2] Observed mean}

\item{cex}{size of points}

\item{HPD}{Highest Posterior Density level for drawing (NA: No HPD drawing )}

\item{xlim}{X-axis range}

\item{ylim}{Y-axis range}

\item{ObsMean}{Draw the observed mean of subgroups if this parameter is TRUE}

\item{...}{other options}
}
\value{
None
}
\description{
plot the posterior response rate with its highest probability density (HPD) interval by subgroup ID
}
\examples{
nDat = c(25, 25, 25, 25) # total number of patients
xDat = c(2, 3, 8, 6)  # number of responses
alpha <- 1e-20 
d0 <- 0.0 
alpha1 = 50   
beta1 = 10  
tau2 <- 0.1  
phi1 <- 0.1  
deltaT <- 0.2  
thetaT <- 0.60   

res <- BCHM(nDat = nDat,
            xDat = xDat,
            alpha = alpha,
            d0 = d0,             
            alpha1 = alpha1, 
            beta1 = beta1,
            tau2 = tau2,
            phi1 = phi1, 
            deltaT = deltaT,
            thetaT = thetaT,
            burnIn = 100,
            MCIter = 200,
            MCNum = 1000,
            seed = 1000
)
print(res$SMatrix)
print(res$Result)
col <- res$Result[,4]

BCHMplot_post_value(res, col, HPD = 0.8)

}
\seealso{
\code{\link{BCHM} Perform the analysis based on the BCHM design.}

\code{\link{BCHMplot_cluster} Plot the clustering results of subgroups. }

\code{\link{BCHMplot_post_dist} Plot the posterior distributions of subgroups. }
}
