diff.plot <- function(diff.before, diff.after, name, weight, title = NULL) {
  # plot standardize difference for before and after weighting
  #
  # Args:
  #   diff.before: standardized mean differece before weighting
  #   diff.after: standardized mean difference after weighting
  #   name: covariate name
  #   weight: propensity score weighting method
  # Returns
  #   A plot is generated

  par(las = 1, lwd = 2, mar = c(5, max(nchar(name)), 4, 2), bty = "n")

  x.range <- range(c(diff.before, diff.after))
  y.range <- c(1, length(name))
  ord <- order(diff.before, decreasing = T)
  plot(
    x = x.range, y = y.range, xaxt = "n", yaxt = "n", type = "n",
    xlim = x.range, ylim = y.range,
    ylab = "", xlab = "Standardized mean difference",
    main = title
  )
  axis(side = 1, at = pretty(x.range))
  axis(side = 2, at = length(name):1, labels = name[ord], tick = F)
  abline(v = 0, lty = 1, col = "gray")
  abline(v = -10, lty = 3, lwd = 1)
  abline(v = 10, lty = 3, lwd = 1)
  points(y = length(name):1, x = diff.before[ord], pch = 4)
  points(y = length(name):1, x = diff.after[ord], pch = 21)
  legend("topleft", legend = c("Unadjusted", weight), pch = c(4, 21))
}

PoolEstimates <- function(Q, U) {
  # Q is point estimate
  # U is variance estimate
  if (length(Q) != length(U)) {
    stop("Q and U are of different length")
  }
  m <- length(Q)
  Q_bar <- mean(Q)
  U_bar <- mean(U) # within-imputation variance
  B <- var(Q) # between-imputation variance
  tt <- U_bar + (1 + 1 / m) * B # total variance
  .tmp <- c(Q_bar, tt)
  names(.tmp) <- c("Q_bar", "total_variance")
  return(.tmp)
}