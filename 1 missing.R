library(tidyverse)

d_path <- file.path(
  "/Users/siyangren/Dropbox (Partners HealthCare)/BOC shared/Chemo during pregnancy (Sella)/data"
)
setwd(d_path)
load(file.path(d_path, "2021-7-19", "data_clean.RData"))

# gestational age at diagnosis missing
df[is.na(gadx), .(chemopg, record_id)]
df %>%
  filter(is.na(gadx)) %>%
  select(record_id, chemopg)
# clean

# breast cancer stage missing
df %>%
  filter(is.na(stage) | stage == 0) %>%
  select(record_id, stage)
# 23, 57, 68, 78, 90

# gestational age at delivery missing
df %>%
  filter(is.na(delga)) %>%
  select(record_id, chemopg)
#   record_id chemopg
# 1        61       1
# 2        63       0
# 3       104       1
# 4       114       0
# 5       121       1
# 6       136       1
# 7       137       1

# preterm birth reason unknown
df %>%
  filter(preterm_birth_reason == "" & ptd == 1) %>%
  select(record_id, sptb, indptb, iatptb, ptd)
# cleans