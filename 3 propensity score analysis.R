packages_list <- c(
  "data.table", "tidyverse", "PSW", "MASS", "sandwich", "mice"
)
lapply(packages_list, library, character.only = TRUE)
select <- dplyr::select

source("/Users/siyangren/Dropbox (Partners HealthCare)/code repository/Functions/R-functions/CreateSummaryTable.R")


d_path <- file.path(
  "/Users/siyangren/Dropbox (Partners HealthCare)/BOC shared/Chemo during pregnancy (Sella)/data"
)
setwd(d_path)
load(file.path("2021-7-19", "data_clean_210915.RData"))


# patients received chemotherapy and gave single birth
df_s <- df %>%
  filter(chemopg == 1 & pregnancy_outcome == "single")

# create a few dummy vars
df_s <- df_s %>%
  mutate(
    age_norm = (agedx - mean(agedx) / sd(agedx)),
    gadx_norm = (gadx - mean(gadx) / sd(gadx)),
    yeardx_c = (yeardx - min(yeardx)) / (max(yeardx) - min(yeardx)),
    I_stage_34 = ifelse(stage > 2, 1, 0), # 12 vs 34
    I_TNBC = ifelse(cancer_type == "TNBC", 1, 0),
    I_HER2 = ifelse(cancer_type == "HER2", 1, 0),
    I_mutation = ifelse(mutation == "Y", 1, 0), # collapse 0 and unknown
    I_neoadjuvant = ifelse(chemo_timing == "neoadjuvant", 1, 0),
    I_surg = surgpg,
    neo_stage = I_neoadjuvant * I_stage_34,
    gadx_surg = gadx_norm * I_surg,
    gadx_stage = gadx_norm * I_stage_34,
    tnbc_surg = I_TNBC * I_surg,
    surg_stage = I_surg * I_stage_34,
    gadx_tnbc = gadx_norm * I_TNBC,
    tnbc_stage = I_TNBC * I_stage_34,
    tnbc_neo = I_TNBC * I_neoadjuvant,
    surg_neo = I_surg * I_neoadjuvant,
    gadx_neo = gadx_norm * I_neoadjuvant
  )

best_f_taxol <- "taxolpreg ~ I_HER2 + gadx_norm + I_TNBC + I_stage_34 +
I_surg + I_neoadjuvant + neo_stage + gadx_tnbc + gadx_surg"

best_f_gcsf <- "gfpreg ~ I_HER2 + gadx_norm + I_TNBC + I_stage_34 +
I_surg + I_neoadjuvant + surg_stage + gadx_tnbc +
tnbc_stage + tnbc_neo"

tbl1_vars <- c(
  "age_norm", "gadx_norm", "I_stage_34",
  "I_TNBC", "I_HER2", "I_neoadjuvant", "I_surg"
)


# Build models
best_model_taxol <- psw(
  data = df_s,
  form.ps = best_f_taxol,
  weight = "MW",
  std.diff = TRUE,
  V.name = tbl1_vars
)

df_s$w_taxol_mw <- best_model_taxol$W

best_model_gcsf <- psw(
  data = df_s,
  form.ps = best_f_gcsf,
  weight = "MW",
  std.diff = TRUE,
  V.name = tbl1_vars
)

df_s$w_gcsf_mw <- best_model_gcsf$W


# model results output
taxol_coefs <- as.data.frame(summary(best_model_taxol$ps.model)$coefficients)
taxol_coefs %>%
  mutate(
    est.or = exp(Estimate),
    lower = exp(Estimate - 1.96 * `Std. Error`),
    upper = exp(Estimate + 1.96 * `Std. Error`),
  ) %>%
  select(est.or, lower, upper, `Pr(>|z|)`)

gcsf_coefs <- as.data.frame(summary(best_model_gcsf$ps.model)$coefficients)
gcsf_coefs %>%
  mutate(
    est.or = exp(Estimate),
    lower = exp(Estimate - 1.96 * `Std. Error`),
    upper = exp(Estimate + 1.96 * `Std. Error`),
  ) %>%
  select(est.or, lower, upper, `Pr(>|z|)`)

# Std diff plots

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


diff.plot(
  diff.before = best_model_taxol$std.diff.before[, "std.diff.pct"],
  diff.after = best_model_taxol$std.diff.after[, "std.diff.pct"],
  name = c(
    "Age",
    "Gestational age",
    "Stage",
    "TNBC",
    "HER2",
    "Neoadjuvant",
    "Surgery"
  ),
  weight = "MW",
  title = "Adjusted for taxane"
)

diff.plot(
  diff.before = best_model_gcsf$std.diff.before[, "std.diff.pct"],
  diff.after = best_model_gcsf$std.diff.after[, "std.diff.pct"],
  name = c(
    "Age",
    "Gestational age",
    "Stage",
    "TNBC",
    "HER2",
    "Neoadjuvant",
    "Surgery"
  ),
  weight = "MW",
  title = "Adjusted for CSF"
)

## Average trt effect models with original dataset, with MI datasets
# update 20210915, combine obstetrical outcome and gestational outcome
# into a single composite outcome
weighted_glm_fit <- function(.df, .outcome, .trt, .weight) {
  .glm <- do.call("glm", list(
    formula = paste0(.outcome, " ~ ", .trt),
    family = binomial,
    data = .df,
    weights = .df[[.weight]]
  ))
  log_odds <- summary(.glm)$coefficients[.trt, "Estimate"]
  log_odds_var <- sandwich(.glm)[.trt, .trt]
  tmp <- c(
    "log_or" = log_odds,
    "var_log_or" = log_odds_var
  )
  return(tmp)
}
# ATE modesl with original dataset
outcome_vars <- c(
  "obstetrical", "gestational", "obstetrical_sens", "gestational_sens"
)
outcome_vars <- c("composite", "composite_sens")
weight_vars <- c("w_taxol_mw", "w_gcsf_mw")
# cbind(df_s$w_taxol_mw, df_s$w_gcsf_mw)
trt_vars <- c("taxolpreg", "gfpreg")

ate_results <- vector("list", length(outcome_vars))
names(ate_results) <- outcome_vars
for (o in outcome_vars) {
  for (i in 1:2) {
    trt <- trt_vars[i]
    w <- weight_vars[i]
    tmp <- weighted_glm_fit(df_s, o, trt, w)
    tmp["est_or"] <- exp(tmp["log_or"])
    tmp["lower_or"] <- exp(tmp["log_or"] - 1.96 * sqrt(tmp["var_log_or"]))
    tmp["upper_or"] <- exp(tmp["log_or"] + 1.96 * sqrt(tmp["var_log_or"]))
    ate_results[[o]][[trt]] <- tmp[c("est_or", "lower_or", "upper_or")]
  }
}

ate_results <- lapply(ate_results, do.call, what = "rbind")

# ATE modesl with MIs
# get multiple imputed datasets first
df_imputed <- sapply(
  outcome_vars,
  function(o) {
    vars <- c(tbl1_vars, o)
    .df_subcols <- df_s %>% select(all_of(vars))
    # mice object
    .fit <- mice(
      data = .df_subcols,
      m = 10, # number of MIs
      maxit = 25, # iteration
      seed = 42,
      method = "pmm"
    )
    # imputed datasets
    imputed <- mice::complete(.fit, "all") # a list
    # add two treatments two the list
    imputed <- lapply(imputed, function(.df) {
      .df$taxolpreg <- df_s$taxolpreg
      .df$gfpreg <- df_s$gfpreg
      return(.df)
    })
    return(imputed)
  },
  simplify = FALSE,
  USE.NAMES = TRUE
)

ate_mi_results <- vector("list", length(outcome_vars))
names(ate_mi_results) <- outcome_vars
# a list, with each outcome as an element
# then each treatment as a sub-element
# each treatment contains a 10*2 matrix, row is iter, col is odds ratio
for (o in outcome_vars) {
  .ls <- df_imputed[[o]]
  for (i in 1:2) {
    trt <- trt_vars[i]
    w <- weight_vars[i]
    placeholder <- matrix(nrow = 10, ncol = 2)
    colnames(placeholder) <- c("log_or", "var_log_or")
    ate_mi_results[[o]][[trt]] <- placeholder
    for (j in 1:10) {
      .df <- .ls[[j]]
      tmp <- weighted_glm_fit(.df, o, trt, w)
      ate_mi_results[[o]][[trt]][j, ] <- tmp
    }
  }
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

# Pool the point estimate and var estimate of all iters
ate_mi_results <- lapply(ate_mi_results, function(.ls) {
  lapply(.ls, function(.matrix) {
    tmp <- PoolEstimates(.matrix[, "log_or"], .matrix[, "var_log_or"])
    tmp["log_or"] <- tmp["Q_bar"]
    tmp["sd_log_or"] <- sqrt(tmp["total_variance"])
    tmp["est_or"] <- exp(tmp["log_or"])
    tmp["lower_or"] <- exp(tmp["log_or"] - 1.96 * tmp["sd_log_or"])
    tmp["upper_or"] <- exp(tmp["log_or"] + 1.96 * tmp["sd_log_or"])
    return(tmp[c("est_or", "lower_or", "upper_or")])
  })
})

# collapse sub-list into a matrix
ate_mi_results <- lapply(ate_mi_results, do.call, what = "rbind")