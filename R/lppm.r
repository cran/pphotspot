#' Fit a Poisson point process model to a point pattern dataset on a linear network
#'
#' Fit a Poisson point process model to a point pattern dataset on a linear network.
#' This function is provided in \pkg{GET} to support \code{\link{hotspots.poislpp}} and
#' \code{\link{hotspots.MatClustlpp}}. See the hotspots vignette available
#' by starting R, typing \code{library("GET")} and \code{vignette("GET")}.
#'
#'
#' The function \code{pois.lppm}, can be used to estimate the inhomogeneous
#' Poisson point process model on linear network. This function provides the
#' \code{firstordermodel}, i.e. the regression model of dependence of crashes on
#' the spatial covariates, \code{EIP}, i.e. estimated inhomogeneous intensity
#' from the data and \code{secondorder}, i.e. estimation of the inhomogeneous
#' K-function. The \code{plot} of the \code{secondorder} provides diagnostics,
#' if the model is adequate for the data. If the estimated $K$-function lies
#' close to the theoretical line, the data does not report any clustering, and
#' the function \code{hotspots.poislpp} can be used for final hotspots detection.
#' If the estimated K-function does not lie close to the theoretical line, and
#' it is above, the data report clustering, and the a clustered point pattern
#' model must be fitted to the data and hotspots detected using this clustered
#' model instead.
#' @param PP Input, a point pattern object (ppp) of spatstat.
#' @param formula An R formula to estimate the first order model.
#'  This formula can contain objects of full size. \code{PP} should be on the
#'  left side of the formula.
#' @param data Data from where the formula takes objects.
#'  Must be acceptable by the function lppm of spatstat.linnet.
#' @param subwin A part of the observation window of \code{PP} to be used for
#' estimating the second order structure. NULL means that the full point pattern
#' is used. Typically this is feasible (not too time consuming).
#' @param r_max The maximum distance on which the K-function is evaluated.
#' Default is computed as \eqn{\sqrt{A}/10}{sqrt(A)/10} where \eqn{A}{A} is the
#' area of the window of observation of \code{X}.
#' @importFrom stats predict
#' @importFrom stats update
#' @importFrom spatstat.linnet lppm
#' @importFrom spatstat.geom area
#' @importFrom spatstat.Knet Knetinhom
#' @return A list of three components:
#' \itemize{
#'   \item EIP = the predicted point process intensity
#'   \item firstordermodel = the model fitted for the intensity
#'   \item secondorder = the inhomogeneous K-function estimated from the given
#' point pattern \code{PP} using the estimated intensity
#' }
#' @export
pois.lppm <- function(PP, formula, data, subwin = NULL, r_max = NULL) {
  PPsub <- PP[, subwin] # Works also for NULL, then PPsub equals PP
  # First order
  M <- lppm(formula, data = data, interaction=NULL)
  EIP <- predict(M)
  # Second order
  A <- spatstat.geom::area(PPsub)
  if(!is.null(r_max)) {
    if(r_max > sqrt(A)/5) warning("Are you sure that your r_max is not too large?")
  }
  else maxr <- sqrt(A)/10
  r <- seq(0, maxr, length.out=100)
  PPsub_K <- Knetinhom(PPsub, lambda = update(M, PPsub), r=r)
  # Return
  list(EIP=EIP, firstordermodel=M, secondorder=PPsub_K)
}

#' Fit a Matern cluster point process to a point pattern dataset on a linear network
#'
#' Fit a Matern cluster point process to a point pattern dataset on a linear network.
#' This function is provided in \pkg{GET} to support
#' \code{\link{hotspots.MatClustlpp}}. See the hotspots vignette available
#' by starting R, typing \code{library("GET")} and \code{vignette("GET")}.
#'
#'
#' The function \code{MatClust.lppm}, can be used to estimate the Matern
#' cluster point pattern with inhomogeneous cluster centers on linear network.
#' This function provides the same outputs as the \code{\link{pois.lppm}} and
#' further estimated parameters \code{alpha} and \code{R}. The \code{secondorder}
#' provides again the diagnostics for checking if the clustered model is
#' appropriate. The sample K-function must be close to the K-function of the
#' estimated model (green line). If it is not the case the searching grid for
#' parameters \code{alpha} and \code{R} that is input in the function must be
#' manipulated to get the a closer result. If the estimated model is adequate
#' one can proceed to the hotspot detection with the use of the function
#' \code{\link{hotspots.MatClustlpp}}. Remark here, that for the estimation of
#' the second order structure a smaller data can be used than for the estimation
#' of the first order structure in order to save the computation time, since the
#' second order is a local characteristics. This smaller window can be specified
#' by \code{subwin}. Then the full point pattern will be used for estimation of
#' first order intensity and the pattern in subwindow will be used for estimating
#' second order characteristic. The input parameters are the same as in
#' \code{\link{pois.lppm}}. Furthermore, \code{valpha}, i.e., vector of proposed
#' alphas which should be considered in the optimization, \code{vR}, i.e., vector
#' of proposed values for R which should be considered in the optimization, must
#' be provided. The user can also specify how many cores should be used in the
#' computation by parameter \code{ncores}.
#' @inheritParams rMatClustlpp
#' @inheritParams pois.lppm
#' @param valpha A vector of parameter values for the parameter alpha of the
#' Matern cluster process.
#' @param vR A vector of parameter values for the parameter R of the Matern
#' cluster process.
#' @param nsim The number of simulated Matern cluster point patterns for
#' evaluating the K-function for any alpha and R values.
#' @param ncores Number of cores used for computations. Default to 1.
#' If NULL, all available cores are used.
#' @importFrom stats update predict
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom parallel clusterEvalQ clusterApplyLB
#' @importFrom spatstat.linnet rpoislpp lppm
#' @importFrom spatstat.Knet Knetinhom
#' @importFrom spatstat.geom area
#' @importFrom GET create_curve_set fdr_envelope
#' @return A list of
#' \itemize{
#'   \item EIP = the predicted point process intensity using the firstordermodel
#'   \item alpha = the estimated value of the parameter alpha of the Matern cluster process, see \code{\link{rMatClustlpp}}
#'   \item R = the estimated value of the parameter R of the Matern cluster process, see \code{\link{rMatClustlpp}}
#'   \item Contrast = the minimum value of the contrast for the estimated alpha and R;
#'   contrast is sqrt(sum((estK - meanK)^2)) where estK is the K-function estimated from the data pattern
#'   and meanK is the mean estimated K function from nsim simulations from the Matern cluster process with
#'   the estimated alpha and R
#'   \item firstordermodel = the model fitted for the intensity using the full
#'   pattern \code{PPfull} (provided to the \code{MatClust.lppm} function in the
#'   argument \code{PP})
#'   \item secondorder = estK, i.e. the K-function estimated from the data pattern
#'   \item MCsecondorder = meanK, i.e. the mean estimated K function from nsim simulations from the Matern cluster process with
#'   the estimated alpha and R
#' }
#' @export
MatClust.lppm <- function(PP, formula, subwin = NULL, valpha, vR, data,
                          nsim = 10, ncores = 1L) {
  PPsub <- PP[, subwin] # Works also for NULL
  lmcppm_par_inner <- function(fakeArg, valpha, EIP, PP, vR, M, r) {
    # Centers from a Poisson process
    Centers <- rpoislpp(EIP/valpha, L=PPsub[['domain']])
    XX <- rMatClustlpp(Centers, vR, valpha, PPsub[['domain']])
    return(Knetinhom(XX, lambda = update(M, XX), r=r)$est)
  }

  if(!is.vector(valpha)) { stop("valpha is not a vector") }
  if(!is.vector(vR)) { stop("vR is not a vector") }
  M <- lppm(formula, data = data)
  EIP <- predict(M)
  maxr <- sqrt(spatstat.geom::area(PP))/10
  r <- seq(0, maxr, length.out=100)
  PPsub_K <- Knetinhom(PPsub, lambda = update(M, PPsub), r=r)
  Contrast <- array(0, c(length(valpha), length(vR)))
  Karray <- array(0, c(length(valpha), length(vR),length(r)))

  fakeArgs <- rep(1, times=nsim)
  # Summarize result of all simulations into one array
  comp_KMC <- function(sim_results, r) {
    KMC <- array(0,length(r))
    for(s in 1:length(sim_results)) {
      KMC <- KMC + sim_results[[s]]
    }
    KMC
  }

  # If the number of cores has not been set by ncores function parameter,
  # try to detect it on the host
  if ( is.null(ncores) ) {
    # detect number of cores
    ncores <- max(1L, detectCores()-1, na.rm = TRUE)
  }

  if(ncores == 1) {
    for(i in 1:length(valpha)) {
      for(j in 1:length(vR)) {
        # Compute the average K of 10 simulation from the model
        KMC_sim_results <- lapply(fakeArgs, lmcppm_par_inner,
                                  valpha[i], EIP, PPsub, vR[j], M, r)
        KMC <- comp_KMC(KMC_sim_results, r)
        # Compute the difference between estimated and average K
        Contrast[i,j] <- sqrt(sum((PPsub_K$est-KMC/nsim)^2))
        Karray[i,j,] <- KMC/nsim
      }
    }
  }
  else {
    # create cluster with ncores nodes
    cl <- makeCluster(ncores, type = "SOCK")
    clusterEvalQ(cl, library("spatstat"))
    clusterEvalQ(cl, library("spatstat.Knet"))
    #print(clusterEvalQ(cl, ls(all.names = TRUE)))

    for(i in 1:length(valpha)) {
      for(j in 1:length(vR)) {
        # Compute the average K of 10 simulation from the model
        KMC_sim_results <- clusterApplyLB(cl, fakeArgs, lmcppm_par_inner,
                                         valpha[i], EIP, PPsub, vR[j], M, r)
        KMC <- comp_KMC(KMC_sim_results, r)
        # Compute the difference between estimated and average K
        Contrast[i,j] <- sqrt(sum((PPsub_K$est-KMC/nsim)^2))
        Karray[i,j,] <- KMC/nsim
      }
    }

    # stop the cluster
    stopCluster(cl)
  }

  # Find the minimum value of the contrast
  id <- which(Contrast == min(Contrast), arr.ind = TRUE)
  alpha <- valpha[id[,1]]
  R <- vR[id[,2]]
  EIP <- EIP/alpha

  list(EIP=EIP,
       alpha=alpha, R=R, Contrast=min(Contrast),
       firstordermodel=M,
       secondorder=PPsub_K, MCsecondorder=Karray[id[,1], id[,2],])
}
