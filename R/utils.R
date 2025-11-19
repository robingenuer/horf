#' Define the regression function used in the simulation by name
#'
#' @param x A vector for one-dimensional functions, a matrix for multidimensional ones, containing values for which the regression function has to be computed
#' @param freg.name Name of the regression function to use, among "sum", "fried1", "sinus" and "abs"
#'
#' @returns A vector of the resulting computed values of the regression function associated to the inputs in `x`.
#' @export
freg <- function(x, freg.name) {
  if (freg.name == "sum") {
    y <- rowSums(x)
  }

  if (freg.name=="fried1") {
    if (length(x[1, ]) < 5) {
      stop("x must have at least 5 coordinates")
    } else {
      y <- 10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 +
        10 * x[, 4] + 5 * x[, 5]
      y <- y / 10
    }
  }

  if (freg.name=="sinus") {
     y <-  sin(2 * pi * x)
  }

  if (freg.name=="abs") {
    y <- abs(x - 0.5)
  }

  return(y)
}


#' Generation of a simulated dataset associated to a given name of the regression function
#'
#' @param nbobs The number of observations of the simulated dataset
#' @param xdim Dimension of the input space of th regression function
#' @param sigma The standard deviation of the centered gaussian noise added to computed output values (default to 1/4)
#' @inheritParams freg
#'
#' @returns
#' @export
#'
#' @examples
#' simData <- simul_data(nbobs = 100, xdim = 1, freg.name = "sinus")
#' curve(freg(x, freg.name = "sinus"), from = 0, to = 1, add = TRUE, col = 2)
simul_data <- function(nbobs, xdim, freg.name, sigma = 1/4) {
  x <- matrix(runif(xdim * nbobs), ncol = xdim)
  eps <- rnorm(nbobs, mean = 0, sd = sigma)
  y <- freg(x, freg.name) + eps
  return(list(x = x, y = y))
}


#' (Internal) Computation of the mean value of the regression function in an hyper-rectangle for the `Fried1` regression function:
#'
#' @param h A matrix with two columns (left and right coordinate of each of the hyper-rectangle dimension) and as many rows as the dimension of the number of parameters of the regression function (at least 5 for the `Fried1` function), that contains the coordinates of an hyper-rectangle associated to a leaf of a tree
#'
#' @returns A vector of length one giving the mean value of the Fried1 regression function in the hyper-rectangle `h`.
#'
beta_fried1 <- function(h) {
  cosc <- function(x) {
    1 / (pi * x) * (cos(pi * h[1, 1] * x) - cos(pi * h[1, 2] * x))
  }
  res <- 10 / ((h[2, 2] - h[2, 1]) * (h[1, 2] - h[1, 1])) *
    stats::integrate(cosc, lower = h[2, 1], upper = h[2, 2])$value +
    20/3 * ((h[3, 2] - 1/2)^2 + (h[3, 2] - 1/2) * (h[3, 1] - 1/2) +
              (h[3, 1] - 1/2)^2) +
    5 * (h[4, 1] + h[4, 2]) + 5/2 * (h[5, 1] + h[5, 2])
  return(res / 10)
}




#' (Internal) Generated a random forest structure build on a simulated dataset using the `randomForest::randomForest()` function
#'
#' @param k Number of leaves of the trees (the same for all trees of a forest)
#' @param seed An optional seed given to the `set.seed()` function
#' @inheritParams simul_data
#' @inheritParams randomForest::randomForest
#'
#' @returns The `forest` component of a `randomForest` class object from the `randomForest` package
#'
#' @examples
#' aforest <- forest_structure(k = 5, xdim = 1, freg.name = "sinus", ntree = 10, nbobs = 100, mtry = 1)
forest_structure <- function(
    k, xdim,  freg.name, ntree, nbobs = NULL, mtry = max(1, floor(d/3)),
    seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  z <- simul_data(nbobs = nbobs, xdim = xdim, freg.name = freg.name)
  x <- z$x
  y <- z$y

  if (xdim == 1) {
    rf <- randomForest::randomForest(
      y ~ x, data = z, ntree = ntree, maxnodes = k, nodesize = 1, mtry = mtry)
  } else {
    rf <- randomForest::randomForest(
      x = x, y = y, ntree = ntree, maxnodes = k, nodesize = 1, mtry = mtry)
  }

  forest <- rf$forest

  return(forest)
}

#' Bounds computation associated to a node of a tree (of a forest)
#'
#' @param forest The forest structure resulting from the `forest_structure` function
#' @param ind_tree Index of the tree to consider
#' @param node Node of the tree to consider
#' @param dim.int a vector of the different dimensions of the input space used for splitting the node, usually (and always for Hold-Out RF) a vector of length one, historically here for dealing for other purely random forests variants
#'
#' @returns A vector of length 2 containing the bounds of the interval associated to that node of that tree
#'
#' @examples
#' aforest <- forest_structure(k = 5, xdim = 1, freg.name = "sinus", ntree = 10, nbobs = 100, mtry = 1)
#' interval_bounds(aforest, 1, 3)
interval_bounds <- function(
    forest, ind_tree, node, dim.int = forest$bestvar[node, ind_tree]) {

  bounds <- matrix(NA, nrow = length(dim.int), ncol = 2)
  left_ancestor <- NULL
  right_ancestor <- NULL
  root <- node

  if (node == 1) {
    bounds[, 1] <- 0
    bounds[, 2] <- 1
  } else {

    while (root != 1) {
      left_ancestor <- c(left_ancestor,
                           which(forest$leftDaughter[, ind_tree] == root))
      right_ancestor <- c(right_ancestor,
                          which(forest$rightDaughter[, ind_tree] == root))
      root <- min(left_ancestor[length(left_ancestor)],
                  right_ancestor[length(right_ancestor)])
    }

    dim_nb <- 1
    for (dim in dim.int) {
      left_ancestor_dim <-
        left_ancestor[which(forest$bestvar[left_ancestor, ind_tree] == dim)]
      right_ancestor_dim <-
        right_ancestor[which(forest$bestvar[right_ancestor, ind_tree] == dim)]

      bsup <- 1
      binf <- 0
      if (length(left_ancestor_dim) > 0) {
        bsup <- forest$xbestsplit[max(left_ancestor_dim), ind_tree]
      }
      if (length(right_ancestor_dim) > 0) {
        binf <- forest$xbestsplit[max(right_ancestor_dim), ind_tree]
      }
      bounds[dim_nb, ] <- c(binf, bsup)
      dim_nb <- dim_nb + 1
    }
  }
  return(bounds)
}


# Hyper-rectangles coordinates computation
hyper_rec_base <- function(ind_tree, forest, hrec, d) {
  hrec.base <- hrec[, ind_tree, , , drop = FALSE]
  nrnodes <- forest$nrnodes
  tree <- lapply(forest[-(1:2)], FUN = function(x)
    subset(x, select = ind_tree, subset = rep(TRUE, nrnodes)))

  for (i in 2: nrnodes) {
    hrec.base[i, , , ] <-
      interval_bounds(
        forest = forest, ind_tree = ind_tree, node = i, dim.int = 1:d)
  }

  return(hrec.base)
}

hyper_rec <- function(forest, d) {
  ntree <- forest$ntree
  nrnodes <- forest$nrnodes
  hrec <- array(data = NA, dim=c(nrnodes, ntree, d, 2))
  hrec[1, , , 1] <- 0
  hrec[1, , , 2] <- 1
  hrec_bases_list <- lapply(
    1:ntree, hyper_rec_base, forest = forest, hrec = hrec, d = d)

  for (j in 1:ntree) {
    hrec[, j, , ] <- hrec_bases_list[[j]]
  }

  return(hrec)
}


# Ideal forest built
forest_tilde <- function(
    k, d, freg.name, ntree, nbobs = NULL, mtry = max(1, floor(d/3)),
    seed = NULL, mc.cores = 1) {

    forest <- forest_structure(
      k = k, d = d, freg.name = freg.name, ntree = ntree,
      nbobs = nbobs, mtry = mtry, seed = seed)

      hrec <- hyper_rec(forest = forest, d = d)

      if (freg.name == "sum") {
        forest$nodepred <- apply(X = hrec, MARGIN = 1:2, FUN = sum) / 2
      }

      if (freg.name == "fried1") {
        forest$nodepred <- apply(X = hrec, MARGIN = 1:2, FUN = beta_fried1)
      }

      if (freg.name == "sinus") {
        forest$nodepred <- matrix(apply(X = hrec, MARGIN = 1:3, FUN = function(u) {
          diff( -cos(2*pi*u) / (2*pi) ) / diff(u)
        }), ncol = ntree)
      }

      if (freg.name == "abs") {
        forest$nodepred <- matrix(apply(X = hrec, MARGIN = 1:3, FUN = function(u) {
          if (u[1] < 0.5 & u[2] > 0.5) {
            out <- (0.5 + u[2]^2 - u[2] + u[1]^2 - u[1]) / (2*u[2] - u[1])
          } else {
            out <- abs((u[1] + u[2] - 1) / 2)
          }
          return(out)
        }), ncol = ntree)
      }

  return(forest)
}


# Descent of an observation in a tree
descent <- function(z, forest, ind_tree,  freg.name,
                    tree_type = "final") {
    node <- 1

    if (tree_type == "final") {
      status <- -1
    }
    if (tree_type == "building") {
      status <- 0
    }

    while (forest$nodestatus[node, ind_tree] != status) {
      if (z[forest$bestvar[node, ind_tree]] <=
          forest$xbestsplit[node, ind_tree]) {
        node <- forest$leftDaughter[node, ind_tree]
      } else {
        node <- forest$rightDaughter[node, ind_tree]
      }
    }
  return(node)
}


# Predict one observation with a forest
f_estim_base <- function(ind_obs, x, forest,  freg.name) {
  pred.trees <- rep(NA, forest$ntree)
  for (ind_tree in 1:forest$ntree) {
    node <- descent(z = x[ind_obs, ], forest = forest, ind_tree = ind_tree,
                    freg.name = freg.name)
    pred.trees[ind_tree] <- forest$nodepred[node, ind_tree]
  }
  return(mean(pred.trees))
}


# Predict all observations of x with a forest
f_estim <- function(x, forest,  freg.name, mc.cores=1) {
    n <- length(x[, 1])
    estims <- unlist(parallel::mclapply(
      1:n, FUN = f_estim_base, x = x, forest = forest,
      freg.name = freg.name, mc.cores=mc.cores))
  return(estims)
}


## TREE
# Ideal tree computation
tree_tilde <- function(
    k, d,  freg.name, nbobs = NULL, mtry = max(1, floor(d/3)),
    seed = NULL) {

    tree <- forest_tilde(
      k = k, d = d,  freg.name = freg.name, ntree = 1,
      nbobs = nbobs, mtry = mtry,
      seed = seed, mc.cores = 1)
  return(tree)
}


# Predict all observations of x with one tree
t_estim <- function(x, tree,  freg.name) {
    estims <- f_estim(
      x = x, forest = tree,  freg.name = freg.name)
  return(estims)
}


# Forest bias computation
bias_for <- function(
    nbleaves,  freg.name, d, nbobs = NULL, beta = 1, nfor = 1,
    r = 500, b = 1, mtry = max(1, floor(d/3)), tree = FALSE,
    var_estim = FALSE, seeds = FALSE,
    mc.cores2 = 1) {

  debut <- Sys.time()
  m <- 1
  for.bias <- vector(mode = "list", length = length(nbleaves))
  for.bias.res <- matrix(NA, nrow = length(nbleaves), ncol = 4)
  x.test <- matrix(runif(d*r), ncol = d)
  s <- freg(x.test, freg.name)
  if (seeds) {seeds <- dget("outputs/seeds")}

  for (k in nbleaves) {
    debutk <- Sys.time()
      x.test.bl <- matrix(runif(d*r, min = b*log(k^(1/d)) / (k^(1/d)),
                                max = 1 - b*log(k^(1/d)) / (k^(1/d))), ncol = d)

    ntree <- ifelse(tree, 1, max(100, floor(k^beta)))

    for.bias[[m]] <- parallel::mclapply(1:nfor, function(j) {
      forest <- forest_tilde(k = k, d = d,
                             freg.name = freg.name, ntree = ntree,
                             nbobs = nbobs, mtry = mtry)
      s.tilde <- f_estim(x = x.test, forest = forest,
                         freg.name = freg.name)
      # mse <- mean((s - s.tilde)^2)
      squaredErrors <- (s - s.tilde)^2

      s.bl <- freg(x.test.bl, freg.name)
      s.tilde.bl <- f_estim(x = x.test.bl, forest = forest,
                            freg.name = freg.name, )
      # mse.bl <- mean((s.bl - s.tilde.bl)^2)
      squaredErrors.bl <- (s.bl - s.tilde.bl)^2

      # return((mse, mse.bl))
      return(cbind(squaredErrors, squaredErrors.bl))
    }, mc.cores = mc.cores2)

    allSquaredErrors <- do.call("rbind", for.bias[[m]])
    for.bias.res[m, c(1, 3)] <- colMeans(allSquaredErrors)

    if (var_estim) {
      allCrossedSquaredErrors <- tcrossprod(allSquaredErrors[, 1])
      for (i in 1:nfor) {
        allCrossedSquaredErrors[1:r + (i-1)*r, 1:r + (i-1)*r] <- NA
      }
      allCrossedSquaredErrors[lower.tri(allCrossedSquaredErrors)] <- NA

      for.bias.res[m, 2] <-
        sqrt(for.bias.res[m, 1]^2 - mean(allCrossedSquaredErrors, na.rm = TRUE))

      allCrossedSquaredErrorsBL <- tcrossprod(allSquaredErrors[, 2])
      for (i in 1:nfor) {
        allCrossedSquaredErrorsBL[1:r + (i-1)*r, 1:r + (i-1)*r] <- NA
      }

      for.bias.res[m, 4] <-
        sqrt(for.bias.res[m, 3]^2 - mean(allCrossedSquaredErrorsBL, na.rm = TRUE))
    }

    # indPairs <- seq(1, nrow(allSquaredErrors), by = r+1)
    # indCrossed <- expand.grid(indPairs, indPairs)
    # allIndCrossed <- indCrossed[-indPairs, ]
    #
    # allIndCrossed <- NULL
    # for (l in 1:10) {
    #   indPairs <- seq(l, nrow(allSquaredErrors), by = r+1)
    #   indCrossed <- expand.grid(indPairs, indPairs)
    #   allIndCrossed <- rbind(allIndCrossed, indCrossed[-indPairs,])
    # }

    # allCrossedSquaredErrors <- t(mapply(function(i, j) {
    #   allSquaredErrors[i,] * allSquaredErrors[j, ]},
    #   allIndCrossed$Var1, allIndCrossed$Var2))
    #
    # for.bias.res[m, c(2, 4)] <-
    #   for.bias.res[m, c(1, 3)]^2 - colMeans(allCrossedSquaredErrors)

    # for.bias.res[m, c(2, 4)] <- apply(for.bias[[m]], MARGIN = 1, sd)
    if (tree) {
      print(paste( "bias tree k =", k, "d =", d, "mtry = ", mtry, "done"))
    } else {
      print(paste( "bias forest k =", k, "d =", d, "mtry = ", mtry, "done"))
    }
    print(Sys.time() - debutk)
    m <- m+1
  }

  names(for.bias) <- nbleaves
  for.bias.res <- cbind(nbleaves, for.bias.res)
  colnames(for.bias.res) <- c("nbleaves", "mean", "sd", "BLmean", "BLsd")
  total_time <- Sys.time() - debut
  time_unit <- attr(total_time, "units")
  print(paste0("Total time = ", total_time))
  save(
    for.bias.res, for.bias, freg.name,  d, nbleaves, nbobs, beta, nfor,
    r, b, mtry, tree, total_time, mc.cores2,
    file = paste0("outputs/",
                  paste(freg.name,  "dim", d,
                        "nbleaves", paste(nbleaves, collapse = "_"),
                        "nbobs", nbobs, "beta", beta, "nfor", nfor, "mtry",
                        mtry, sep="_"),
                  ifelse(tree, "_TREE", "_FOREST"),
                  "_BIAS.Rdata"))
  return(for.bias)
}

# Slope computation
slope <- function(y, x) {
  data <- data.frame(x, y)
  return(lm(y ~ x, data = data)$coef[2])
}

# Intercept computation
intercept <- function(y, x) {
  data <- data.frame(x, y)
  return(lm(y ~ x, data = data)$coef[1])
}


CIslope <- function(
    for.bias, nbleaves,  freg.name, d, nbobs = NULL, beta = 1,
    tree = FALSE, nfor = 1, nbSamples = 500, alpha = 0.05) {

  resamplings <- lapply(1:nbSamples, function(i) {
    indU <- sample(
      1:length(for.bias[[1]]), size = length(for.bias[[1]]), replace = TRUE)
    indX <- sample(1:nrow(for.bias[[1]][[1]]), size = nrow(for.bias[[1]][[1]]),
                   replace = TRUE)
    bias <- t(sapply(for.bias, function(biask) {
      subBiask <- lapply(biask[indU], function(subB) { subB[indX, ] })
      subSquaredErrorsk <- do.call("rbind", subBiask)
      return(colMeans(subSquaredErrorsk))
    }))
    log.bias <- log(bias, base = 2)
    log.nbleaves <- log(nbleaves, base = 2)
    slopes <- apply(log.bias, MARGIN = 2, FUN = slope, x = log.nbleaves)
    intercepts <- apply(log.bias, MARGIN = 2, FUN = intercept, x = log.nbleaves)
    return(list(bchaps = bias, slopes = slopes, intercepts = intercepts))
  })

  resampledBchaps <- sapply(resamplings, function(resamp) { resamp$bchaps[, 1] })
  bchapsCI <- t(apply(resampledBchaps, MARGIN = 1, FUN = function(bchapsk) {
    sortedBchapsk <- sort(bchapsk)
    bchapsBounds <- sortedBchapsk[c(floor(nbSamples * alpha / 2),
                                    ceiling(nbSamples * (1 - alpha / 2)))]
    return(bchapsBounds)
  }))

  resampledBchapsBL <- sapply(resamplings, function(resamp) { resamp$bchaps[, 2] })
  bchapsBLCI <- t(apply(resampledBchapsBL, MARGIN = 1, FUN = function(bchapsk) {
    sortedBchapsk <- sort(bchapsk)
    bchapsBounds <- sortedBchapsk[c(floor(nbSamples * alpha / 2),
                                    ceiling(nbSamples * (1 - alpha / 2)))]
    return(bchapsBounds)
  }))

  resampledSlopes <- sapply(resamplings, function(resamp) { resamp$slopes })
  slopesCI <- t(apply(resampledSlopes, MARGIN = 1, FUN = function(slopes) {
    sortedSlopes <- sort(slopes)
    bounds <- sortedSlopes[c(floor(nbSamples * alpha / 2),
                             ceiling(nbSamples * (1 - alpha / 2)))]
    return(bounds)
  }))

  resampledIntercepts <- sapply(resamplings, function(resamp) { resamp$intercepts })
  interceptsCI <- t(apply(resampledIntercepts, MARGIN = 1, FUN = function(intercepts) {
    sortedIntercepts <- sort(intercepts)
    bounds <- sortedIntercepts[c(floor(nbSamples * alpha / 2),
                                 ceiling(nbSamples * (1 - alpha / 2)))]
    return(bounds)
  }))
  load(file = paste0("outputs/",
                     paste(freg.name,  "dim", d,
                           "nbleaves", paste(nbleaves, collapse = "_"),
                           "nbobs", nbobs, "beta", beta, "nfor", nfor, "mtry",
                           mtry, sep="_"),
                     ifelse(tree, "_TREE", "_FOREST"),
                     "_BIAS.Rdata"))
  save(for.bias.res, for.bias, freg.name,  d, nbleaves, nbobs, beta,
       nfor, r, b, mtry, tree, total_time, mc.cores2, slopesCI,
       interceptsCI, bchapsCI, bchapsBLCI,
       file = paste0("outputs/",
                     paste(freg.name,  "dim", d,
                           "nbleaves", paste(nbleaves, collapse = "_"),
                           "nbobs", nbobs, "beta", beta, "nfor", nfor, "mtry",
                           mtry, sep="_"),
                     ifelse(tree, "_TREE", "_FOREST"),
                     "_BIAS.Rdata"))
  return(list(slopesCI = slopesCI, interceptsCI = interceptsCI,
              bchapsCI = bchapsCI, bchapsBLCI = bchapsBLCI))
}

