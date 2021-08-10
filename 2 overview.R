library(tidyverse)

source("C:/Users/sren/Dropbox (Partners HealthCare)/code repository/Functions/R-functions/CreateSummaryTable.R")

d_path <- file.path(
  "C:/Users/sren/Dropbox (Partners HealthCare)/BOC shared/Chemo during pregnancy (Sella)/data"
)
setwd(d_path)

load(file.path("2021-7-19", "data_clean.RData"))

# --- factor vars need to be factor class --- #
# table 1 factor vars
tbl_1_factor_vars <- c(
  "dxtri", "stage", "cancer_type", "mutation", "chemo_timing",
  "chemo", "gcsf", "surgpg"
)
df <- df %>%
  mutate(across(all_of(tbl_1_factor_vars), as.factor))
# table 2 factor vars
tbl_2_factor_vars <- c(
  "pregnancy_outcome", "preterm_birth_outcome", "preterm_birth_reason"
)
df <- df %>%
  mutate(across(all_of(tbl_2_factor_vars), as.factor))
# table 3 factor vars
tbl_3_factor_vars <- c(
  "sptb", "pprom", "gdm", "pih", "chorio",
  "sga", "apgar5v7", "anomalies"
)
df <- df %>%
  mutate(across(all_of(tbl_3_factor_vars), as.factor))


# --- Generate tables --- #

# Table 1
tbl1_vars <- c(
  "agedx", "gadx", "year_cat", "dxtri", "stage", "cancer_type",
  "mutation", "chemoga", "chemo_timing", "chemo", "gcsf", "surgpg"
)
summaryTable(df, tbl1_vars, "chemopg", .digits = 1)
summaryTable(df %>% filter(chemopg == 1), tbl1_vars, "taxolpreg", .digits = 1)


# Table 2
tbl2_vars <- c(
  "delga", "pregnancy_outcome",
  "preterm_birth_outcome", "preterm_birth_reason"
)
summaryTable(df, tbl2_vars, "chemopg", .digits = 1)

# Table 3
# all twin+ pregnancies should be excluded from table 3, 4
tbl3_vars <- c(
  "obstetrical_sens", "sptb", "pprom", "gdm", "pih", "chorio",
  "gestational_sens", "sga", "apgar5v7", "anomalies"
)

# single birth + chemopg == 1
df_sub <- df %>% filter(chemopg == 1 & pregnancy_outcome == "single")
# table 3, split by whether they got taxol during pregnancy
summaryTable(df_sub, tbl3_vars, "taxolpreg", .digits = 1)

# Table 4
tbl4_vars <- c(
  "obstetrical_sens", "sptb", "pprom", "gdm", "pih", "chorio",
  "gestational_sens", "sga", "apgar5v7", "anomalies"
)
summaryTable(df_sub, tbl4_vars, "gfpreg", .digits = 1)