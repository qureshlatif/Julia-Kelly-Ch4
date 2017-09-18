library(dplyr)

setwd("F:/research stuff/FS_PostDoc/consult_&_collaborate/JuliaH_dissert/Chap4")  # Set workspace to location with 'Data_compiled' workspace.
load("Data_compiled.RData")

tab.sum <- Tree.data %>% tbl_df %>%
  filter(site_name %in% c("tv_treat", "tv_cont") & species == "pien" & year == 2015 |
           site_name %in% c("cannibal", "slumpass") & species == "pien" & year == 2014) %>%
  filter(status %in% c("healthy", "early Inf")) %>%
  filter(!is.na(dbh)) %>%
  group_by(status) %>%
  summarise(Mean = mean(dbh), SD = sd(dbh), Min = min(dbh), Max = max(dbh), n = n())

