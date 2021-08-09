packages_list <- c(
  "data.table", "tidyverse", "kableExtra", "PSW", "MASS", "sandwich", "mice"
)
lapply(packages_list, library, character.only = TRUE)
select <- dplyr::select

source("/Users/siyangren/Dropbox (Partners HealthCare)/code repository/Functions/R-functions/CreateSummaryTable.R")


d_path <- file.path(
  "/Users/siyangren/Dropbox (Partners HealthCare)/BOC shared/Chemo during pregnancy (Sella)/data"
)
setwd(d_path)

load(file.path("2021-7-19", "data_clean.RData"))

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

## MI
mice_fit_ls <- sapply(
  c("obstetrical", "gestational", "obstetrical_sens", "gestational_sens"),
  function(outc) {
    vars <- c(tbl1_vars, outc)
    .sub_data <- df_s %>% select(all_of(vars))
    # mice object
    .fit <- mice(
      data = .sub_data,
      m = 10, # number of MIs
      maxit = 25, # iteration
      seed = 42,
      method = "pmm"
    )
    # imputed datasets
    imputed <- mice::complete(.fit, "all") # a list
    # add two treatments two the list
    imputed <- lapply(imputed, function(.data) {
      .data$taxolpreg <- df_s$taxolpreg
      .data$gfpreg <- df_s$gfpreg
      return(.data)
    })
    return(imputed)
  },
  simplify = FALSE,
  USE.NAMES = TRUE
)

# Average treatment effect model
ATE_composite <- sapply(
  c("obstetrical", "gestational", "obstetrical_sens", "gestational_sens"),
  function(outc) {
    sapply(c("taxolpreg", "gfpreg"), function(trt) {
      w <<- if (trt == "taxolpreg") df_s$w_taxol_mw else df_s$w_gcsf_mw
      # outc ~ trt
      .glm <- glm(
        formula = paste0(outc, " ~ ", trt),
        family = binomial,
        data = df_s,
        weights = w
      )
      # extract coefs
      .log_odds <- summary(.glm)$coefficients[trt, "Estimate"]
      .log_odds_var <- sandwich(.glm)[trt, trt]
      return(data.frame(
        est.log.or = .log_odds,
        std.log.or = sqrt(.log_odds_var)
      ))
    }, simplify = FALSE, USE.NAMES = TRUE) %>%
      rbindlist(., idcol = "Treatment")
  },
  simplify = FALSE, USE.NAMES = TRUE
) %>%
  rbindlist(., idcol = "Outcome") %>%
  mutate(
    est.or = exp(est.log.or),
    lower = exp(est.log.or - 1.96 * std.log.or),
    upper = exp(est.log.or + 1.96 * std.log.or)
  ) %>%
  arrange(Treatment, Outcome)

ATE_composite_MI <- sapply(
  c("obstetrical", "gestational", "obstetrical_sens", "gestational_sens"),
  function(outc) {
    .ls <- mice_fit_ls[[outc]]
    # tmp is a dataframe, with 10 iterations and 2 treatments
    tmp <- lapply(.ls, function(.data) {
      # taxane
      tmp1 <- {
        .glm <- glm(
          paste0(outc, " ~ taxolpreg"),
          family = binomial,
          data = .data,
          weights = df_s$w_taxol_mw
        )
        .log_odds <- summary(.glm)$coefficients["taxolpreg", "Estimate"]
        .log_odds_var <- sandwich(.glm)["taxolpreg", "taxolpreg"]
        # return
        data.frame(
          Treatment = "taxolpreg",
          log.or = .log_odds,
          log.or.var = .log_odds_var
        )
      }
      # CSF
      tmp2 <- {
        .glm <- glm(
          paste0(outc, " ~ gfpreg"),
          family = binomial,
          data = .data,
          weights = df_s$w_gcsf_mw
        )
        .log_odds <- summary(.glm)$coefficients["gfpreg", "Estimate"]
        .log_odds_var <- sandwich(.glm)["gfpreg", "gfpreg"]
        # return
        data.frame(
          Treatment = "gfpreg",
          log.or = .log_odds,
          log.or.var = .log_odds_var
        )
      }
      # combine taxane and CSF
      return(
        rbind(tmp1, tmp2)
      )
    }) %>%
      rbindlist(., idcol = "Iter")
    # combine stats over 10 imputed datasets
    tmp <- tmp %>%
      group_by(Treatment) %>%
      summarise(
        est.log.or = PoolEstimates(log.or, log.or.var)[1],
        std.log.or = sqrt(PoolEstimates(log.or, log.or.var)[2])
      ) %>%
      mutate(
        est.or = exp(est.log.or), # convert from coef to odds ratio
        lower = exp(est.log.or - 1.96 * std.log.or),
        upper = exp(est.log.or + 1.96 * std.log.or)
      ) %>%
      as.data.frame()
    return(tmp)
  },
  simplify = FALSE,
  USE.NAMES = TRUE
) %>%
  rbindlist(., idcol = "Outcome") %>%
  mutate(Outcome = paste0(Outcome, "_mi")) %>%
  arrange(Outcome, Treatment)

ATE_overall <- rbind(
  ATE_composite %>%
    select(Outcome, Treatment, est.or, lower, upper),
  ATE_composite_MI %>%
    select(Outcome, Treatment, est.or, lower, upper)
)

# For taxane
# taxane and obstetrical
ATE_overall %>%
  filter(Treatment == "taxolpreg" & str_detect(Outcome, "^obs")) %>%
  select(Outcome, est.or, lower, upper)
# taxane and gestational
ATE_overall %>%
  filter(Treatment == "taxolpreg" & str_detect(Outcome, "^ges")) %>%
  select(Outcome, est.or, lower, upper)
# GCSF and obstetrical
ATE_overall %>%
  filter(Treatment == "gfpreg" & str_detect(Outcome, "^obs")) %>%
  select(Outcome, est.or, lower, upper)
# GCSF and gestational
ATE_overall %>%
  filter(Treatment == "gfpreg" & str_detect(Outcome, "^ges")) %>%
  select(Outcome, est.or, lower, upper)