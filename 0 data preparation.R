library(tidyverse)

d_path <- file.path(
  "C:/Users/sren/Dropbox (Partners HealthCare)/BOC shared/Chemo during pregnancy (Sella)/data"
)
setwd(d_path)
df <- read.csv(
  file.path("2021-7-19", "RegistryOfPregnancyA_DATA_2021-07-19_1531.csv")
)

# exclude record ID 26, 56, 128
names(df)[1] <- "record_id"
df <- df %>%
  filter(!record_id %in% c(26, 56, 128))
# gestational age at checmotherapy initiation is not numeric
df <- df %>% mutate(
  chemoga = as.numeric(case_when(
    chemoga %in% c("n/a", "unknown") ~ as.character(NA),
    chemoga == "18 weeks" ~ "18",
    TRUE ~ chemoga
  ))
)

# table 1 changes
# Impute gestational age at diagnosis by gestational age at surgery and
# weeks between diagnosis and surgery
df <- df %>%
  mutate(gadx = ifelse(
    is.na(gadx),
    as.numeric(surgga) - as.numeric(time2surg),
    gadx
  )) %>%
  mutate(gadx = ifelse(gadx < 0, 0, gadx))
# Categorize year of diagnosis
df <- df %>%
  mutate(
    year_cat = case_when(
      yeardx < 2000 ~ "before 2000",
      2000 <= yeardx & yeardx < 2005 ~ "2000-2005",
      2005 <= yeardx & yeardx < 2010 ~ "2005-2010",
      2010 <= yeardx & yeardx < 2015 ~ "2010-2015",
      2015 <= yeardx & yeardx < 2020 ~ "2015-2020",
      yeardx >= 2020 ~ "2020+"
    ),
    year_cat = factor(
      year_cat,
      levels = c(
        "before 2000", "2000-2005", "2005-2010",
        "2010-2015", "2015-2020", "2020+"
      ),
      ordered = TRUE
    )
  )
# Stage 0 is missing
df <- df %>%
  mutate(stage = ifelse(stage == 0, NA, stage))
# cancer type
df <- df %>%
  mutate(cancer_type = case_when(
    her2 == 1 ~ "HER2",
    her2 == 0 & erpr == 1 ~ "HR",
    her2 == 0 & erpr == 0 ~ "TNBC"
  ))
# breast cancer gene
df <- df %>%
  mutate(mutation = case_when(
    brca %in% c("BRCA 1", "BRCA1", "BRCA2", "PABL2") ~ "Y",
    brca %in% c("0", "no", "No") ~ "N",
    brca %in% c("", "not tested") ~ "Unknown"
  ))
# chemo timing (could be out of pregnancy)
df <- df %>%
  mutate(chemo_timing = case_when(
    anychemo == 1 & neoadjchemo == 1 ~ "neoadjuvant",
    anychemo == 1 & (neoadjchemo == 0 | is.na(neoadjchemo)) ~ "adjuvant",
    anychemo == 0 ~ "missing"
  ))
# GCSF
df <- df %>%
  mutate(gcsf = case_when(
    gfpreg == 0 ~ "no",
    gfpreg == 1 & gftype == "daily GCSF" ~ "daily_gcsf",
    TRUE ~ "depot_gcsf"
  ))
# Chemotherapy
df <- df %>%
  mutate(chemo = case_when(
    chemopg == 0 ~ "no chemo",
    acpreg == 1 & taxolpreg == 0 ~ "ac_only",
    acpreg == 1 & taxolpreg == 1 ~ "ac_taxane",
    acpreg == 0 & taxolpreg == 1 ~ "taxane_only"
  ))
# dose densed AC
df <- df %>% rename(dose_dense_ac = ddacpreg)


# table 2 changes
# pregnancy outcome
df <- df %>%
  mutate(
    pregnancy_outcome = case_when(
      twin == 0 ~ "single birth",
      twin == 1 ~ "double birth",
      outcomepg == "TAB" ~ "therapeutic abortion",
      outcomepg == "SAB" ~ "spontaneous abortion",
      TRUE ~ as.character(NA)
    ),
    preterm_birth = case_when(
      pregnancy_outcome %in% c("single birth", "double birth") &
        ptd == 1 ~ "Y",
      pregnancy_outcome %in% c("single birth", "double birth") &
        ptd == 0 ~ "N",
      pregnancy_outcome %in% c("single birth", "double birth") &
        is.na(ptd) ~ "unknown",
      !pregnancy_outcome %in% c("single birth", "double birth") ~ "NA"
    ),
    preterm_reason = case_when(
      preterm_birth == "Y" & sptb == 1 ~ "spontaneous",
      preterm_birth == "Y" & indptb == 1 ~ "Medically indicated",
      preterm_birth == "Y" & iatptb == 1 ~ "Iatrogenic",
      TRUE ~ "NA"
    )
  )



# Table 3 changes
df <- df %>%
  mutate(
    apgar5 = as.numeric(ifelse(apgar5 == "6; 7", "6", apgar5)),
    apgar5v7 = as.numeric(apgar5 < 7)
  )



# Composite outcomes
# update 20210915: add a single composite outcome, exclude gestational
# hypertesion and gestational diabetes
df$composite <- apply(
  df %>% select(sptb, pprom, chorio, sga, apgar5v7, anomalies),
  1,
  function(x) {
    ifelse(all(is.na(x)), NA, as.numeric(sum(x, na.rm = TRUE) > 0))
  }
)
df$composite_sens <- apply(
  df %>% select(sptb, pprom, chorio, sga, apgar5v7, anomalies),
  1,
  function(x) {
    ifelse(sum(x, na.rm = TRUE) > 0, # if any individual outcome is 1
      1, # define as 1
      ifelse(sum(is.na(x)) > 0, # else if any is NA
        NA, # define as NA
        0
      )
    )
  }
)
# four options
# primary analysis: composite outcome is missing only if
# all individual outcomes are missing

df$obstetrical <- apply(
  df %>% select(sptb, pprom, gdm, pih, chorio),
  1,
  function(x) {
    ifelse(all(is.na(x)), # if all invidual outcomes are NA
      NA, # define as NA
      as.numeric(sum(x, na.rm = TRUE) > 0)
    ) # else will be either 0 or 1
  }
)
# sensitivity analysis: composite outcome is missing if any
# individual outcome is missing
df$obstetrical_sens <- apply(
  df %>% select(sptb, pprom, gdm, pih, chorio),
  1,
  function(x) {
    ifelse(sum(x, na.rm = TRUE) > 0, # if any individual outcome is 1
      1, # define as 1
      ifelse(sum(is.na(x)) > 0, # else if any is NA
        NA, # define as NA
        0
      )
    )
  }
)
df$gestational <- apply(
  df %>% select(sga, apgar5v7, anomalies),
  1,
  function(x) {
    ifelse(all(is.na(x)),
      NA,
      as.numeric(sum(x, na.rm = TRUE) > 0)
    )
  }
)
df$gestational_sens <- apply(
  df %>% select(sga, apgar5v7, anomalies),
  1,
  function(x) {
    ifelse(sum(x, na.rm = TRUE) > 0,
      1,
      ifelse(sum(is.na(x)) > 0,
        NA,
        0
      )
    )
  }
)

# check the outcomes
df %>% distinct(sptb, pprom, gdm, pih, chorio, obstetrical, obstetrical_sens)

df %>%
  distinct(
    sptb, pprom, chorio, sga, apgar5v7, anomalies, composite, composite_sens
  )

save(df, file = file.path(d_path, "2021-7-19", "data_clean_211022.RData"))