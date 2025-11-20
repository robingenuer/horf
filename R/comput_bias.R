#' Forest bias computation
#'
#' @param nbleaves A vector of the different number of leaves considered to compute the forest bias at
#' @param nfor Number of forests to build (default to 1)
#' @param nbobs_test Number of observations of the simulated test set (default to 500)
#' @param tree A logical value: TRUE for a tree bias computation, FALSE for a forest bias computation (defaut to FALSE)
#' @param var_estim A logical value to decide to estimate the variance of bias computed for the `nfor` forests built
#' @param mc.cores Number of cores to give to `parallel::mclapply()` function
#' @inheritParams forest_tilde
#'
#' @returns A list containing all bias computations. The first component is a matrix summarizing those results with mean and variance estimations for each `nbleaves` value. The second argument is a list containing all individual results and the last component is the call of the function.
#' @export
#' @importFrom parallel mclapply
#'
#' @examples
#' bias_result <- comput_bias(nbleaves = c(2^5, 2^6), freg.name = "sinus", xdim = 1, nbobs = 640, nbobs_test = 100, nfor = 2, var_estim = TRUE)
comput_bias <- function(
    nbleaves,  freg.name, xdim, nbobs = NULL, nfor = 1, nbobs_test = 500,
    mtry = max(1, floor(xdim/3)), tree = FALSE, var_estim = FALSE,
    seeds = FALSE, mc.cores = 1) {

  debut <- Sys.time()
  cl <- match.call()
  cl[[1]] <- as.name("comput_bias")
  m <- 1
  for.bias <- vector(mode = "list", length = length(nbleaves))
  for.bias.res <- matrix(NA, nrow = length(nbleaves), ncol = 4)
  x.test <- matrix(runif(xdim*nbobs_test), ncol = xdim)
  s <- freg(x.test, freg.name)
  if (seeds) {seeds <- dget("outputs/seeds")}

  for (k in nbleaves) {
    debutk <- Sys.time()
    x.test.bl <- matrix(
      runif(xdim*nbobs_test, min = log(k^(1/xdim)) / (k^(1/xdim)),
            max = 1 - log(k^(1/xdim)) / (k^(1/xdim))),
      ncol = xdim)

    ntree <- ifelse(tree, 1, max(100, k))

    for.bias[[m]] <- parallel::mclapply(1:nfor, function(j) {
      forest <- forest_tilde(k = k, xdim = xdim,
                             freg.name = freg.name, ntree = ntree,
                             nbobs = nbobs, mtry = mtry)
      s.tilde <- f_estim(zmat = x.test, forest = forest)
      squaredErrors <- (s - s.tilde)^2

      s.bl <- freg(x.test.bl, freg.name)
      s.tilde.bl <- f_estim(zmat = x.test.bl, forest = forest)
      squaredErrors.bl <- (s.bl - s.tilde.bl)^2

      return(cbind(squaredErrors, squaredErrors.bl))
    }, mc.cores = mc.cores)

    allSquaredErrors <- do.call("rbind", for.bias[[m]])
    for.bias.res[m, c(1, 3)] <- colMeans(allSquaredErrors)

    if (var_estim) {
      allCrossedSquaredErrors <- tcrossprod(allSquaredErrors[, 1])
      for (i in 1:nfor) {
        allCrossedSquaredErrors[1:nbobs_test + (i-1)*nbobs_test, 1:nbobs_test + (i-1)*nbobs_test] <- NA
      }
      allCrossedSquaredErrors[lower.tri(allCrossedSquaredErrors)] <- NA

      for.bias.res[m, 2] <-
        sqrt(for.bias.res[m, 1]^2 - mean(allCrossedSquaredErrors, na.rm = TRUE))

      allCrossedSquaredErrorsBL <- tcrossprod(allSquaredErrors[, 2])
      for (i in 1:nfor) {
        allCrossedSquaredErrorsBL[1:nbobs_test + (i-1)*nbobs_test, 1:nbobs_test + (i-1)*nbobs_test] <- NA
      }

      for.bias.res[m, 4] <-
        sqrt(for.bias.res[m, 3]^2 - mean(allCrossedSquaredErrorsBL, na.rm = TRUE))
    }

    if (tree) {
      print(paste( "bias tree k =", k, "xdim =", xdim, "mtry = ", mtry, "done"))
    } else {
      print(paste( "bias forest k =", k, "xdim =", xdim, "mtry = ", mtry, "done"))
    }
    print(Sys.time() - debutk)
    m <- m+1
  }

  names(for.bias) <- nbleaves
  for.bias.res <- cbind(nbleaves, for.bias.res)
  colnames(for.bias.res) <- c("nbleaves", "mean", "sd", "BLmean", "BLsd")
  total_time <- Sys.time() - debut
  time_unit <- attr(total_time, "units")
  print(paste("Total time =", total_time, time_unit))
  res <- list(for.bias.res = for.bias.res, for.bias = for.bias, call = cl)

  return(res)
}
