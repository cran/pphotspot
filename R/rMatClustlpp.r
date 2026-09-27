#' Simulating the Matern cluster process on a linear network
#'
#' Simulating the Matern clusters with parameters alpha and R given the (x,y)-
#' coordinates of the parent points.
#' This function is provided in \pkg{GET} to support
#' \code{\link{hotspots.MatClustlpp}}. See the hotspots vignette available
#' by starting R, typing \code{library("GET")} and \code{vignette("GET")}.
#'
#' @param Centers The (x,y)-coordinates of parent points
#' @param R Cluster radius parameter of the Matern cluster process.
#' @param alpha Parameter related to mean number of points per cluster.
#' @param LL The linear network on which the point pattern should be simulated.
#' @param check_vol Logical. TRUE for checking if the ball producing the cluster
#'  has any intersection with linear network.
#' @importFrom sf st_linestring st_sfc st_buffer st_intersection
#' @importFrom spatstat.geom psp Window lengths_psp disc volume
#' @importFrom spatstat.random rpoisppOnLines
#' @importFrom spatstat.linnet lpp
#' @return A point pattern on linear network, \code{lpp} object of spatstat.
#' @examples
#' # example code
#' if(require(spatstat.geom, quietly=TRUE) & require(spatstat.linnet, quietly=TRUE)) {
#' data("roadcrash")
#' win <- owin(xrange = roadcrash$xrange,
#'             yrange = roadcrash$yrange)
#' X <- ppp(x = roadcrash$x, y = roadcrash$y, window = win)
#' Vertices.pp <- ppp(x = roadcrash$Vertices.x,
#'                    y = roadcrash$Vertices.y,
#'                    window=win)
#' L <- linnet(vertices=Vertices.pp,
#'             edges = roadcrash$Edges)
#' PPfull <- lpp(X, L)
#' subwin <- owin(c(-760000, -740000), c(-1160000, -1140000))
#' PPsub <- PPfull[, subwin]
#' Centers <- rpoislpp(0.0005, L=PPsub[['domain']])
#' XX <- rMatClustlpp(Centers, 500, 5, PPsub[['domain']])
#' plot(Centers, col=1, pch=16)
#' points(XX, col=2)
#' }
#' @export
rMatClustlpp <- function(Centers, R, alpha, LL, check_vol=FALSE) {
  e <- as.matrix(LL$lines$ends)
  e2 <- list()
  for(i in 1:nrow(e)) {
    e2[[i]] <- st_linestring(matrix(e[i,], 2,2, byrow=TRUE))
  }
  LN <- sf::st_sfc(e2)
  X <- array(0,0)
  Y <- array(0,0)
  cs <- st_sfc(apply(cbind(Centers$data$x, Centers$data$y), 1,
                         sf::st_point, simplify=FALSE))
  bufs <- st_buffer(cs, R, nQuadSegs = 8)
  LN1s <- st_intersection(bufs, sf::st_union(LN))
  stopifnot(length(LN1s)==length(Centers$data$x))
  for (p in seq_along(Centers$data$x)) {
    segment_matrix <- geometry_to_segments(LN1s[p])

    if (nrow(segment_matrix) == 0L) {
      next
    }

    LL1 <- psp(
      segment_matrix[, "x0"],
      segment_matrix[, "y0"],
      segment_matrix[, "x1"],
      segment_matrix[, "y1"],
      window = Window(LL),
      check = FALSE
    )

    vol <- sum(lengths_psp(LL1))
    if(check_vol) {
      BBCOutD_ss <- disc(radius=R, centre=c(Centers$data$x[p],
                                            Centers$data$y[p]),
                         npoly = 32)
      vol2 <- volume(LL[BBCOutD_ss])
      if(abs(vol-vol2) > 1e-4) stop("vol")
    }
    Xp <- rpoisppOnLines(alpha/vol, L=LL1)
    X <- append(X, as.numeric(Xp$x))
    Y <- append(Y, as.numeric(Xp$y))
  }
  lpp(cbind(X,Y), LL)
}
