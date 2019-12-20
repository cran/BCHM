
#' @title Perform the analysis based on the BCHM design. 
#'
#' @description The BCHM computation is based on the Bayesian Cluster Hierarchical Model (BCHM) to perform the non-parametric Bayesian clustering and posterior probability calculation with the Bayesian hierarchical model for binary response data in several subgroups. Due to the potential heterogeneity among subgroups, the exchangeability assumption across all subgroups may not hold. A Bayesian non-parametric method is applied to calculate the number of clusters by conducting the multiple cluster classification based on subgroup outcomes. Due to the MCMC sampling, the clustering result is dynamic.  A cluster matrix (Similarity Matrix) is constructed to depict the strength of association between any two subgroups to be classified into the same cluster. The Bayesian hierarchical model is used to compute the posterior probability of treatment effect with the borrowing strength determined by the similarity matrix values calculated from the Bayesian non-parametric clustering.
#' @param nDat Number of patients in each subgroup.
#' @param xDat Number of responses in each subgroup.
#' @param mu Hyperprior mean for the cluster. 
#' @param sigma02 Hyperprior variance for the cluster.
#' @param sigmaD2 Variance of subgroup response rate.
#' @param alpha Alpha value of the Dirichlet Process determining number of clusters.
#' @param d0  Minimum value for the similarity matrix.
#' @param alpha1 Prior for borrowing strength gamma(alpha1, beta1) in the hierarchical model.
#' @param beta1 Prior for borrowing strength gamma(alpha1, beta1) in the hierarchical model.
#' @param tau2 Hyperprior precision parameter of subgroup means in the hierarchical model
#' @param phi1 The response rate of the standard treatment.
#' @param deltaT The expected improvement in the response rate over the standard treatment.
#' @param thetaT Threshold value for the inference claiming efficacy.
#' @param burnIn Number of burn_in in MCMC.
#' @param MCIter Number of MCMC iterations.
#' @param MCNum Number of MCMC iterations in the hierarchical model.
#' @param seed  Random seed.
#' @return The return is a list including three elements: \code{Samples}, \code{SMatrix}, and \code{Result}.
#' @return The return list element \code{Samples} is the sampled posterior response rates of all subgroups.
#' @return The return list element \code{SMatrix} is the similarity matrix of all subgroups.
#' @return The return list element \code{Result} is the calculation results of all subgroups. It has seven columns: the number of responses of each subgroups, the number of patients in each subgroups, the observe response rates of each subgroups, the cluster index of each subgroups, the posterior mean response rates of each subgroups, the probability of Pr(P>Phi1+deltaT) of each subgroups, and the final decision (0: Not rejected the null, 1 Rejected the null). Note: Because a MCMC computation is applied in the clustering calculation, there are many possible clustering configurations. The cluster index in column 4 is the most possible clustering configuration. 
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
#' BCHMplot_post_value(res, col, HPD = 0.8)
#' BCHMplot_post_dist(res, col, lty=1:length(nDat), lwd =3, xlim=c(0, 0.8))
#' 
#' @seealso \code{\link{BCHMplot_cluster} Plot the clustering results of subgroups. }
#' @seealso \code{\link{BCHMplot_post_value} Plot the posterior response of subgroups. }
#' @seealso \code{\link{BCHMplot_post_dist} Plot the posterior distributions of subgroups. }
#' @importFrom stats median var sd runif density rnorm
#' @importFrom graphics plot
#' @importFrom rjags jags.model coda.samples dic.samples
#' 
#' @export
BCHM <- function(nDat,
                 xDat,
                 mu = 0.2, 
                 sigma02 = 10, 
                 sigmaD2 = 0.001,                 
                 alpha = 1e-60,
                 d0 = 0.05,
                 alpha1 = 50,
                 beta1 = 10,
                 tau2 = 0.1,
                 phi1 = 0.1,
                 deltaT = 0.05,
                 thetaT = 0.6,                  
                 burnIn = 10000,
                 MCIter = 20000,
                 MCNum = 20000,
                 seed = 1000)
{
  numArm <- length(nDat)

  if (numArm != length(xDat))
  {
    stop("Numbers of subgroups in nDat and xDat are not equal.")
  }
  if (numArm > 20)
  {
    stop("Numbers of subgroups is more than 20.")
  }
  
  set.seed(seed)
  weight <- nDat
  delta <- deltaT
  alphaP <- alpha
  alpha <- alpha1
  beta <- beta1
  
  posi <- xDat
  phi2 <- phi1 + delta
  rate <- c()

  priorMean <- mu        # Prior of the mean
  priorVar <- sigma02    # Prior of the variance
  
  res <- posi / weight
  x <- as.matrix(res)
  estVarGroup <- sigmaD2   # Data variance
  inSD <- sd(res)
  #print("Computation starts running:")
  result <-
    gibbsSampler(x,
                 alphaP,
                 priorMean,
                 priorVar,
                 estVarGroup,
                 weight,
                 burnIn,
                 MCIter)
  tables <- result$tables
  
  sm <- matrix(0, numArm, numArm)
  tSize <- dim(tables)[1]
  for (i in 1:tSize)
  {
    rr <- tables[i,]
    for (j in 1:numArm)
    {
      for (k in 1:numArm)
      {
        if (rr[j] == rr[k])
        {
          sm[j, k] <- sm[j, k] + 1
        }
      }
    }
  }
  
  
  sm <- sm / tSize
  smMinR <- d0
  sm[sm < smMinR] <- smMinR
  
  #browser()
  
  #print("The similarity Matrix C_ij:")
  
  SMatrix <- sm
  rownames(SMatrix) <- 1:numArm
  colnames(SMatrix) <- 1:numArm
  #print(SMatrix)
  

  optm <- optimizeSil(x, tables)
  ind <- optm$ind
  #perfom<-optm$sil
  #ind <-1
  perform <- 0.0

  clusters <- result$tables[ind,]
  nCluster <- max(clusters)
  outAll <- matrix(0, numArm, 7)
  colnames(outAll) <-
    c("Resp.",
      "No.Pat",
      "Obs. Rate",
      "Cluster",
      "Post. Mean",
      "Prob(P>Phi1+deltaT)",
      "Decision") #, "Independent Mean")

  outAll[, 1] <- posi
  outAll[, 2] <- weight
  outAll[, 3] <- round(posi / weight + 1e-9, 3)
  outAll[, 4] <- clusters
  #outAll[,7]<-rate
  
  #  browser()
  sm[sm < 0.001] <- 0.001 # Prevent bugs crash
  maxY <- 0
  #print("Computing the posterior distribution of all subgroups:")

  for (i in 1:numArm)
  {

    cInd <- 1:numArm
    t <- length(cInd)
    smP <- sm[i,]
    
    mydata <-
      list(
        y = posi[cInd],
        n = weight[cInd],
        m = smP,
        numGroups = length(cInd),
        targetResp = phi2,
        mu0 = logit(mean(posi[cInd] / weight[cInd])),
        tau2 = tau2,
        alpha = alpha,
        beta = beta
      )
    
    
    mText1 <- modelStr()
    modelSpec1 <-textConnection(mText1)
    
    parameters <- c("p")
    
    jSeed <- floor(runif(1, 1, 10000))
    
    jags1 <- jags.model(
      modelSpec1,
      data = mydata,
      n.chains = 4,
      n.adapt = MCNum / 4,
      quiet = TRUE,
      inits = list(.RNG.name = "base::Wichmann-Hill",
                   .RNG.seed = jSeed)
    )
    
    MCRes1 <- coda.samples(jags1,
                           parameters,
                           n.iter = MCNum,
                           verbose = FALSE,
                           progress.bar = "none",
                           thin = 1)

    samples <- MCRes1[[1]]  

    sampledP <- samples[, i] 

    maxY <- max(maxY, max(density(sampledP)$y))
    if (i == 1)
    {
      allPost <- sampledP
    } else{
      allPost <- cbind(allPost, sampledP)
    }

    index <- cInd[i]
    prob <- sum(sampledP > phi2) / length(sampledP)
    outAll[index, 6] <- round(prob, 3)
    outAll[index, 5] <- round(mean(sampledP), 3)
  }
  
  allDat <-
    list(
      outAll = outAll,
      clusterFit = perform,
      sm = sm,
      clusterT = tables
    )
  #dput(allDat,"result.txt")
  #allDat
  
  decision <- outAll[,6] > thetaT
  outAll[,7] <- decision
  rownames(outAll) <- 1:numArm
  colnames(allPost) <- 1:numArm

  for(ii in 1:numArm)
  {
    rownames(outAll)[ii] <- paste0("Subg. ", ii)
  }  
  result <- list(Samples = allPost, SMatrix = SMatrix, Result = outAll)
  class(result) <-"BCHM_result"
  result
  
  #outAll[,6]
}