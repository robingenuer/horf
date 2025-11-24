library(ggplot2)
library(scales)


#' Title
#'
#'
#' @param forest.comput.bias.res Resulting object from `comput_bias()` function for forests (with `tree = FALSE` in the call of the function)
#' @param forest.confint.res Resulting object from `comput_confint()` function applied to the content of `forest.comput.bias.res`
#' @param tree.comput.bias.res Resulting object from `comput_bias()` function for trees (with `tree = TRUE` in the call of the function)
#' @param tree.confint.res Resulting object from `comput_confint()` function applied to the content of `tree.comput.bias.res`
#'
#' @returns Nothing is returned, only a plot is made.
#' @export
#'
#' @import ggplot2
#'
#' @examples
#' forest_bias_result <- comput_bias(nbleaves = c(2^5, 2^6), freg.name = "sinus", xdim = 1, nbobs = 640, nbobs_test = 100, nfor = 2, var_estim = TRUE)
#' forest_confint_result <- comput_confint(forest_bias_result, nbBootstrap = 50)
#' tree_bias_result <- comput_bias(nbleaves = c(2^5, 2^6), freg.name = "sinus", xdim = 1, nbobs = 640, nbobs_test = 100, nfor = 2, tree = TRUE, var_estim = TRUE)
#' tree_confint_result <- comput_confint(tree_bias_result, nbBootstrap = 50)
#' plot_bias(forest_bias_result, tree_bias_result)
plot_bias <- function(forest.comput.bias.res, forest.confint.res, tree.comput.bias.res, tree.confint.res) {

  for.bias <- data.frame(forest.comput.bias.res$for.bias.res)
  for.slopesCI <- forest.confint.res$slopesCI
  for.interceptsCI <- forest.confint.res$interceptsCI
  for.bchapsCI <- forest.confint.res$biasesCI

  tree.bias <- data.frame(tree.comput.bias.res$for.bias.res)
  tree.slopesCI <- tree.confint.res$slopesCI
  tree.interceptsCI <- tree.confint.res$interceptsCI
  tree.bchapsCI <- tree.confint.res$biasesCI

  if (!all.equal(tree.bias[, 1], for.bias[, 1])) {
    stop("number of leaves must be the same for tree and forest")
  } else {
    bias <- cbind(tree.bias[, -1], for.bias[, -1])
  }
  means <- seq(1, ncol(bias), 2)
  sds <- seq(2, ncol(bias), 2)
  log.bias <- log(bias, base = 2)
  log.nbleaves <- log(tree.bias$nbleaves, base = 2)

  pentes <- apply(log.bias[, means], MARGIN = 2,
                  FUN = function(y, x = log.nbleaves) {
                    return(lm(y ~ x, data = data.frame(x, y))$coef[2])})
  pentes <- round(pentes, 3)

  intercepts <- apply(log.bias[, means], MARGIN = 2,
                      FUN = function(y,x = log.nbleaves) {
                        return(lm(y ~ x, data = data.frame(x, y))$coef[1])})
  # consts <- round(2^intercepts, 3)
  intercepts <- round(intercepts, 3)

  slopesCI <- round(rbind(tree.slopesCI, for.slopesCI), 3)
  rownames(slopesCI) <- c("treeSquaredErr", "forestSquaredErr")

  interceptsCI <- round(rbind(tree.interceptsCI, for.interceptsCI), 3)
  rownames(interceptsCI) <- c("treeSquaredErr", "forestSquaredErr")

  tree.bias[, 3:4] <- tree.bchapsCI
  colnames(tree.bias)[3:4] <- c("Binf", "Bsup")
  tree.bias$Type <- "tree"
  for.bias[, 3:4] <- for.bchapsCI
  colnames(for.bias)[3:4] <- c("Binf", "Bsup")
  for.bias$Type <- "forest"
  bias <- rbind(tree.bias, for.bias)

  pbiasBars <- ggplot(data = bias, aes(x = nbleaves, y = mean, color = Type)) +
    geom_point() + geom_line() +
    geom_errorbar(aes(ymin = Binf, ymax = Bsup), width = 0.2) +
    scale_y_continuous(trans = log2_trans(),
                       breaks = trans_breaks("log2", function(x) 2^x),
                       labels = trans_format("log2", math_format(2^.x))) +
    scale_x_continuous(trans = log2_trans(),
                       breaks = trans_breaks("log2", function(x) 2^x),
                       labels = trans_format("log2", math_format(2^.x))) +
    scale_color_manual(
      values = viridis::viridis(3)[1:2],
      breaks = c("tree", "forest"),
      labels = c(paste0("tree", "\n",
                        "r=", pentes[1], "(",
                        paste(slopesCI[1, ], collapse = ","), ")", "\n",
                        "i=", intercepts[1], "(",
                        paste(interceptsCI[1, ], collapse = ","), ")"),
                 paste0("forest", "\n",
                        "r=", pentes[2], "(",
                        paste(slopesCI[2, ], collapse = ","), ")" , "\n",
                        "i=", intercepts[2], "(",
                        paste(interceptsCI[2, ], collapse = ","), ")"))) +
    xlab("Number of leaves of trees (log-scale)") + ylab("Bias (log-scale)") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

  print(pbiasBars)
}

