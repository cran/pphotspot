# Hot spots detection from Poisson or Matern Cluster point pattern
# - parallel version
#' @importFrom stats predict
#' @importFrom spatstat.linnet rpoislpp lppm density.lpp
#' @importFrom GET fdr_envelope
hotspots_par <- function(model=c("pois", "MatClust"),
                         PP, formula, data,
                         clusterparam=NULL,
                         sigma = 250, nsim = 10000,
                         ncores = 1L, ...) {
  if(model == "pois") {
    hotspots_par_inner <- function(fakeArg, clusterparam, EIP, PP, sigma) {
      simss <- rpoislpp(EIP, L=PP[['domain']])
      return(density.lpp(simss, sigma = sigma, distance="euclidean"))
    }
  }
  else { # "MatClust", clusterparam must be provided
    hotspots_par_inner <- function(fakeArg, clusterparam, EIP, PP, sigma) {
      Centers <- rpoislpp(EIP, L=PP[['domain']])
      simss <- rMatClustlpp(Centers, clusterparam[['R']],
                            clusterparam[['alpha']], PP[['domain']])
      return(density.lpp(simss, sigma = sigma, distance="euclidean"))
    }
  }

  M <- lppm(formula, data = data)
  if(model == "pois") EIP <- predict(M)
  else EIP <- predict(M)/clusterparam[['alpha']]
  densi <- density.lpp(PP, sigma = sigma, distance="euclidean")
  sims.densi <- vector(mode = "list", length = nsim)
  fakeArgs <- rep(1, times=nsim)

  # if the number of cores has not been set by ncores function parameter,
  # try to detect it on the host
  if(is.null(ncores)) {
    # detect number of cores
    ncores <- max(1L, detectCores()-1, na.rm = TRUE)
  }

  if(ncores == 1) {
    sims.densi <- lapply(fakeArgs, hotspots_par_inner,
                         clusterparam, EIP, PP, sigma)
  }
  else {
    # create cluster with ncores nodes
    cl <- makeCluster(ncores, type = "SOCK")
    clusterEvalQ(cl, library("spatstat"))
    sims.densi <- clusterApplyLB(cl, fakeArgs, hotspots_par_inner,
                                 clusterparam, EIP, PP, sigma)
    # stop the cluster
    stopCluster(cl)
  }

  yx <- expand.grid(densi$yrow, densi$xcol)

  noNA_id <- which(!is.na(densi$v))
  noNA_idsim <- which(!is.na(sims.densi[[1]]$v))
  noNA_id <- intersect(noNA_id, noNA_idsim)
  #max(densi[noNA_id]); summary(sapply(sims.densi, FUN=max))

  cset <- create_curve_set(list(
    r=data.frame(x=yx[,2], y=yx[,1],
                 width=densi$xstep, height=densi$ystep)[noNA_id,],
    obs=as.vector(densi$v)[noNA_id],
    sim_m=sapply(sims.densi, FUN=function(x){ as.vector(x$v)[noNA_id] },
                 simplify=TRUE)))
  res <- fdr_envelope(cset, alternative = "greater", ...)
  res
}

#' Hotspot detection for Poisson point process
#' - parallel version
#'
#' See the hotspots vignette available by starting R, typing
#' \code{library("GET")} and \code{vignette("GET")}.
#' @param PP The point pattern living in a network.
#' @param formula A formula for the intensity.
#' @param sigma To be passed to density.lpp.
#' @param nsim Number of simulations to be performed.
#' @param ... Additional parameters to be passed to \code{\link[GET]{fdr_envelope}}.
#' @inheritParams MatClust.lppm
#' @references Mrkvička et al. (2023). Hotspot detection on a linear network in the presence of covariates: A case study on road crash data. DOI: 10.2139/ssrn.4627591
#' @return An object of type \code{fdr_envelope} and \code{global_envelope},
#' see further details in \code{\link[GET]{global_envelope_test}}.
#' @export
hotspots.poislpp <- function(PP, formula, data,
                             sigma=250, nsim = 10000,
                             ncores=1L, ...) {
  hotspots_par(model="pois", PP=PP, formula=formula, data=data,
               subwin=NULL, # No need for subwin, fast
               sigma=sigma, nsim=nsim, ncores=ncores, ...)
}

#' Hotspot detection for Matern point process
#' - parallel version
#'
#' See the hotspots vignette available by starting R, typing
#' \code{library("GET")} and \code{vignette("GET")}.
#' @inheritParams pois.lppm
#' @inheritParams hotspots.poislpp
#' @inheritParams rMatClustlpp
#' @references Mrkvička et al. (2023). Hotspot detection on a linear network in the presence of covariates: A case study on road crash data. DOI: 10.2139/ssrn.4627591
#' @return An object of type \code{fdr_envelope} and \code{global_envelope},
#' see further details in \code{\link[GET]{global_envelope_test}}.
#' @export
hotspots.MatClustlpp <- function(PP, formula, R, alpha, data,
                                 sigma = 250, nsim = 10000,
                                 ncores = 1L, ...) {
  hotspots_par(model="MatClust", PP=PP, formula=formula, data=data,
               clusterparam=list(alpha=alpha, R=R),
               sigma=sigma, nsim=nsim, ncores=ncores, ...)
}
