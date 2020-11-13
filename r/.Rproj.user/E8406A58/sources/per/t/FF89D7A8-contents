if (! "data.table" %in% installed.packages()[, 'Package']) {
  install.packages('data.table')
}

library(data.table)

data_path <- "../data/2020-8-17/"
redcap_data <- fread(paste0(data_path, "RegistryOfPregnancyA_DATA_2020-08-17_0951.csv"))
names(redcap_data)[1] <- "record_id"
redcap_data <- redcap_data %>% filter(! record_id %in% c(26, 56, 128)) # exclude record ID 26, 56, 128

# gestational age at checmotherapy initiation is not numeric
redcap_data$chemoga[redcap_data$chemoga == "n/a"] <- NA
redcap_data$chemoga[redcap_data$chemoga == "unknown"] <- NA
redcap_data$chemoga[redcap_data$chemoga == "18 weeks"] <- "18"
redcap_data$chemoga <- as.numeric(redcap_data$chemoga)

# table 1 changes
# cancer type
redcap_data[her2 == 1, cancer_type := "HER2"]
redcap_data[her2 == 0 & erpr == 1, cancer_type := "HR"]
redcap_data[her2 == 0 & erpr == 0, cancer_type := "TNBC"]
# breast cancer gene
redcap_data[brca %in% c("BRCA 1", "BRCA1", "BRCA2", "PABL2"), mutation := "Y"]
redcap_data[brca %in% c("0", "no", "No"), mutation := "N"]
redcap_data[brca %in% c("", "not tested"), mutation := "Unknown"]
# chemo timing (could be out of pregnancy)
redcap_data[anychemo == 1 & neoadjchemo == 1, chemo_timing := "neoadjuvant"]
redcap_data[anychemo == 1 & neoadjchemo == 0, chemo_timing := "adjuvant"]
# GCSF
redcap_data[gfpreg == 0, gcsf := "no"]
redcap_data[gfpreg == 1 & gftype == "daily GCSF", gcsf := "daily_gcsf"]
redcap_data[gfpreg == 1 & gftype %in% c("depot GCSF", "Neulasta", "Pegfilgrastim"), gcsf := "depot_gcsf"]
# Chemotherapy
redcap_data[chemopg == 0, chemotherapy := "no chemo"]
redcap_data[acpreg == 1 & taxolpreg == 0, chemotherapy := "AC_only"]
redcap_data[acpreg == 1 & taxolpreg == 1, chemotherapy := "AC_taxane"]
redcap_data[acpreg == 0 & taxolpreg == 1, chemotherapy := "taxane_only"]

# table 2 changes
# pregnancy outcome
redcap_data[twin == 0, pregnancy_outcome := "single"]
redcap_data[twin == 1, pregnancy_outcome := "double"]
redcap_data[outcomepg == "TAB", pregnancy_outcome := "Termination"]
redcap_data[outcomepg == "SAB", pregnancy_outcome := "Spontaneous"]
redcap_data[is.na(twin) & is.na(outcomepg), pregnancy_outcome := NA]
# preterm birth
redcap_data[ptd == 0, preterm_birth_outcome := ">=37"]
redcap_data[ptd == 1, preterm_birth_outcome := "<37"]
redcap_data[sptb == 1, preterm_birth_reason := "Spontaneous"]
redcap_data[indptb == 1, preterm_birth_reason := "Medical"]
redcap_data[sptb == 1 & indptb == 1, preterm_birth_reason := "Spontaneous & Medical"]
redcap_data[iatptb == 1, preterm_birth_reason := "Iatrogenic"]
redcap_data[sptb == 0 & indptb == 0 & iatptb == 0, preterm_birth_reason := "Unknown"]
redcap_data[ptd == 0, preterm_birth_reason := "Not applicable"]

# Table 3 changes
redcap_data$apgar5[redcap_data$apgar5 == "6; 7"] <- "6" # ?
redcap_data$apgar5 <- as.numeric(redcap_data$apgar5)
redcap_data[, apgar5v7 := as.numeric(apgar5 < 7)]
# table 3 composite outcomes
redcap_data[, obstetrical := as.numeric(rowSums(.SD, na.rm = TRUE) > 0), .SDcols = c("sptb", "pprom", "gdm", "pih", "chorio")]
redcap_data[, gestational := as.numeric(rowSums(.SD, na.rm = TRUE) > 0), .SDcols = c("sga", "apgar5v7", "anomalies")]

write.csv(redcap_data, file = paste0(data_path, "cleaned_data.csv"), row.names = FALSE)

