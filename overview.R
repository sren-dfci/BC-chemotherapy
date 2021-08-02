library(tidyverse)

source("/Users/siyangren/Dropbox (Partners HealthCare)/code repository/Functions/R-functions/CreateSummaryTable.R")

d_path <- file.path(
  "/Users/siyangren/Dropbox (Partners HealthCare)/BOC shared/Chemo during pregnancy (Sella)/data"
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
  mutate(across(tbl_1_factor_vars, as.factor))
# table 2 factor vars
tbl_2_factor_vars <- c(
  "pregnancy_outcome", "preterm_birth_outcome", "preterm_birth_reason"
)
df <- df %>%
  mutate(across(tbl_2_factor_vars, as.factor))
# table 3 factor vars
tbl_3_factor_vars <- c(
  "sptb", "pprom", "gdm", "pih", "chorio",
  "sga", "apgar5v7", "anomalies"
)
df <- df %>%
  mutate(across(tbl_3_factor_vars, as.factor))


# --- Generate tables --- # 
df_chemo <- df %>% filter(chemopg == 1)
df_no_chemo <- df %>% filter(chemopg == 0)
df_taxane <- df %>% filter(taxolpreg == 1)
df_no_taxane <- df %>% filter(taxolpreg == 0)


