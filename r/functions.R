numStats <- function(x, label = NULL, .digits = 2, useNA = "always") {
  # returns a 2 columns dataframe

  if ("data.frame" %in% class(x)) {
    if (ncol(x) > 1) stop("too many columns in x")
    label <- names(x)
    x <- x[[label]] # change x from df to vector
  }
  useNA <- match.arg(useNA)
  # median, 25, 75 quantile
  funs <- c(mean, sd)
  m_sd <- round(sapply(funs, function(f) f(x, na.rm = TRUE)), .digits)
  x_miss <- is.na(x)
  .miss <- round(c(sum(x_miss), mean(x_miss) * 100), .digits)
  if (useNA == "always" | (useNA == "ifany" & .miss[1] > 0)) {
    .output <- data.frame(
      Variable = c("mean (SD)", "missing"),
      Value = c(
        paste0(m_sd[1], " (", m_sd[2], ")"),
        paste0(.miss[1], " (", .miss[2], "%)")
      )
    )
  } else {
    .output <- data.frame(
      Variable = "mean (SD)",
      Value = paste0(m_sd[1], " (", m_sd[2], ")")
    )
  }
  if (!is.null(label)) {
    .output <- rbind(
      data.frame(Variable = paste0(label, "  "), Value = ""),
      .output
    )
  }
  return(.output)
}

factorStats <- function(x, label = NULL, .digits = 2, useNA = "always") {
  if ("data.frame" %in% class(x)) {
    if (ncol(x) > 1) stop("too many columns in x")
    label <- names(x)
    x <- x[[label]]
  }
  useNA <- match.arg(useNA)
  if (useNA == "ifany") {
    x <- addNA(x, ifany = TRUE)
  } else if (useNA == "always") {
    x <- addNA(x, ifany = FALSE)
  }
  tbl <- as.data.frame(table(x))
  tbl_prop <- as.data.frame(round(prop.table(table(x)) * 100, .digits),
    responseName = "Percentage"
  )
  merged_tbl <- merge(tbl, tbl_prop, by = names(tbl)[1], all = TRUE)
  names(merged_tbl)[1] <- "Variable"
  merged_tbl$Value <- paste0(merged_tbl$Freq, " (", merged_tbl$Percentage, "%)")
  merged_tbl <- merged_tbl[order(merged_tbl$Variable), c("Variable", "Value")]
  if (!is.null(label)) {
    merged_tbl <- rbind(
      data.frame(Variable = paste0(label, "  "), Value = ""),
      merged_tbl
    )
  }
  return(merged_tbl)
}

calibration.plot <- function(.obs, .pred, .nbreaks = 10, .breaks = NULL, .type = c("equal.length", "quantiles"), .title = "Calibration Plot", .xlim = NULL, .ylim = NULL, .points = F,
                             .xlab = "Predicted probability", .ylab = "Fraction with actual outcome", ...) {
  # Calibration plot
  # Plot observed proportions vs. predicted probabilities, binned by predicted probability groups
  # Visualization of Hosmer-Lemeshow test
  # See Clinical Prediction Models, E. Steyerberg, p271
  # Cut predictions into L quantiles or L regions of equal probability
  # Plot midpoints vs. proportion of cases with outcome in that group

  # rm(.obs, .pred, .nbreaks, .type, .title, .xlim, .ylim, .points, .xlab, .ylab, BREAKS, .pred.CUT, Y, X)

  if (is.null(.xlim)) .xlim <- range(.pred)
  if (is.null(.ylim)) .ylim <- range(.pred)
  plot(0, 0, col = "white", xlim = .xlim, ylim = .ylim, xlab = .xlab, ylab = .ylab)

  if (.points) points(jitter(.pred, 0.05), jitter(.obs, 0.05), pch = 1)

  if (.type[1] == "quantiles") {
    print(paste("Using", .nbreaks, "equal n intervals (quantiles) for breaks"))
    BREAKS <- quantile(.pred, probs = seq(0, 1, l = .nbreaks))
  } else if (.type[1] == "equal.length") {
    print(paste("Using", .nbreaks, "equal x-axis length intervals for breaks"))
    BREAKS <- seq(min(.pred), max(.pred), l = .nbreaks)
  } else if (!is.null(.breaks)) {
    print("Using user-specified bins")
    BREAKS <- .breaks
  }
  .pred.CUT <- cut(.pred, breaks = BREAKS, include.lowest = T, ordered_result = T)
  # table(.pred.CUT, exclude = NULL)

  if (is.factor(.obs)) .obs <- as.numeric(.obs == "1") # Assumes the event is 1

  Y <- tapply(.obs, list(.pred.CUT), mean)
  N <- tapply(.obs, list(.pred.CUT), length)
  SE <- tapply(.obs, list(.pred.CUT), sd) / sqrt(N)
  if (.type[1] == "quantiles") {
    X <- tapply(.pred, list(.pred.CUT), median)
    # This plots at the median value - more appropriate for quantile-based breaks...
    # But is this cheating? I do not think so...
    # You are comparing the observed prevalence in a bin to the median phat in that bin.
    # For a rare outcome, the right-most bin distribution of phat is quite skewed, so comparing to the mid-point or mean is unfair.
  } else {
    X <- as.numeric(BREAKS[-length(BREAKS)] + diff(BREAKS) / 2) # This plots at the midpoint
  }

  points(X, Y, pch = 2)
  for (i in 1:length(X)) {
    ERROR <- Y[i] + c(-1, 1) * 1.96 * SE[i]
    ERROR[ERROR < 0] <- 0
    ERROR[ERROR > 1] <- 1
    lines(rep(X[i], 2), ERROR)
  }
  abline(a = 0, b = 1, lty = 2)
  SUB <- !is.na(X) & !is.na(Y)
  lines(smooth.spline(X[SUB], Y[SUB], df = 3))
  title(.title)


  lm_model <- lm(Y ~ X, data = data.frame(Y = Y, X = X))

  return(data.frame(
    intercept = round(lm_model$coefficients[1], 3),
    slope = round(lm_model$coefficients[2], 3)
  ))
}

prob.hist <- function(.obs, .pred, .xlim = NULL, .breaks = seq(0, 1, 0.01), .cutoff = NULL, ...) {
  if (is.null(.xlim)) .xlim <- range(.pred)
  par(mfrow = c(3, 1))
  hist(.pred, breaks = .breaks, xlim = .xlim, ...)
  if (!is.null(.cutoff)) abline(v = .cutoff, col = "green")
  abline(v = median(.pred), col = "red")
  hist(.pred[.obs == 0], breaks = .breaks, xlim = .xlim, ...)
  abline(v = median(.pred[.obs == 0]), col = "red")
  if (!is.null(.cutoff)) abline(v = .cutoff, col = "green")
  hist(.pred[.obs == 1], breaks = .breaks, xlim = .xlim, ...)
  abline(v = median(.pred[.obs == 1]), col = "red")
  if (!is.null(.cutoff)) abline(v = .cutoff, col = "green")
  par(mfrow = c(1, 1))
  return("Histogram of predicted probabilities")
}

StdDiff <- function(.grp, .var, type = "numeric") {
  # https://support.sas.com/resources/papers/proceedings12/335-2012.pdf
  # .grp should be 0/1
  # .var shoule be either numeric of factor
  # type could be any of numeric, binary and category

  if (length(unique(.grp)) != 2) stop(".grp should have two levels")
  if (length(.grp) != length(.var)) stop(".grp and .var have different length")

  .grp <- as.factor(.grp)
  if (type %in% c("binary", "category")) .var <- droplevels(as.factor(.var)) # drop unused levels

  xc <- .var[.grp == 0] # control group
  xt <- .var[.grp == 1] # trt group
  if (sum(is.na(xc))) {
    xc <- na.omit(xc)
    print(sum(is.na(xc)), " missing in control group")
  }
  if (sum(is.na(xt))) {
    xt <- na.omit(xt)
    print(sum(is.na(xt)), " missing in treatment group")
  }
  # d
  if (type == "numeric") {
    .d <- abs(mean(xc) - mean(xt)) / sqrt((sd(xc)^2 + sd(xt)^2) / 2) # (m1-m2) / sqrt((var1 + var2) / 2)
  } else if (type == "binary") {
    P0 <- mean(xc == levels(xc)[1])
    P1 <- mean(xt == levels(xt)[1])
    .d <- abs(P0 - P1) / sqrt((P0 * (1 - P0) + P1 * (1 - P1)) / 2)
  } else if (type == "category") {
    P0 <- as.matrix(as.data.frame(prop.table(table(xc)))$Freq)[-1, , drop = FALSE] # N-1 * 1 matrix
    P1 <- as.matrix(as.data.frame(prop.table(table(xt)))$Freq)[-1, , drop = FALSE]
    S <- matrix(rep(0, (nlevels(.var) - 1)^2), ncol = nlevels(.var) - 1)
    for (i in 1:(nlevels(.var) - 1)) {
      for (j in 1:(nlevels(.var) - 1)) {
        if (i == j) {
          S[i, j] <- (P0[i] * (1 - P0[i]) + P1[i] * (1 - P1[i])) * 0.5
        } else {
          S[i, j] <- (P0[i] * P0[j] + P1[i] * P1[j]) * -0.5
        }
      }
    }
    .d <- sqrt(t(P0 - P1) %*% solve(S) %*% (P0 - P1))[1, 1]
  }
  # sd(d)
  nc <- length(xc)
  nt <- length(xt)
  se <- sqrt((nc + nt) / (nc * nt) + .d^2 / (nc + nt) * 0.5)
  .result <- c(.d, se, .d - 1.96 * se, .d + 1.96 * se)
  names(.result) <- c("d", "sd(d)", "lower CI", "upper CI")
  .result <- round(.result, 3)
  # output
  # print(paste0(type, " standardized difference"))
  return(.result)
}

extractATE <- function(.psw) {
  .wt_estimator <- .psw[
    c(
      "est.risk.wt", "std.risk.wt", "est.rr.wt", "std.rr.wt",
      "est.or.wt", "std.or.wt"
    )
  ] # list
  .wt_estimator <- data.frame(
    matrix(
      unlist(.wt_estimator),
      ncol = length(.wt_estimator)
    )
  )
  colnames(.wt_estimator) <- c(
    "est.risk.wt", "std.risk.wt", "est.rr.wt", "std.rr.wt",
    "est.or.wt", "std.or.wt"
  )
  .wt_estimator
}

extractStdDiff <- function(.psw) {
  .df <- data.frame(
    Vars = row.names(.psw$std.diff.before),
    Before = unname(.psw$std.diff.before[, 5]),
    After = unname(.psw$std.diff.after[, 5])
  )
  .df
}

# https://github.com/cran/PSW/blob/master/R/psw-internal.R
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

## MI estimates of treatment
# pool MI results
# https://stats.idre.ucla.edu/wp-content/uploads/2016/02/multipleimputation.pdf
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