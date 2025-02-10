pacman::p_load(here, tidyverse, data.table)

source(file.path("/ihome/spurkayastha/soumik/2025_biometrika/code/analysis_element/2025_01_26_bpd.R"))
source(file.path("/ihome/spurkayastha/soumik/2025_biometrika/code/analysis_element/2025_01_26_bps.R"))

dia <- read_csv(file.path("/ihome/spurkayastha/soumik/2025_biometrika/output/element_2/bpd_summary.csv")) %>% 
  mutate(across(contains("_"), ~1000*(as.numeric(.x)))) %>% 
  mutate(path1 = paste0(sprintf('%0.3f', fit_g1), " (", sprintf('%0.3f', lcb_g1), ", ", sprintf('%0.3f', ucb_g1), ")"), 
         path2 = paste0(sprintf('%0.3f', fit_g2), " (", sprintf('%0.3f', lcb_g2), ", ", sprintf('%0.3f', ucb_g2), ")")) %>% 
  select(type, g1, g2, path1, path2) %>% 
  arrange(type, g1, g2)


sys <- read_csv(file.path("/ihome/spurkayastha/soumik/2025_biometrika/output/element_2/bps_summary.csv")) %>% 
  mutate(across(contains("_"), ~1000*(as.numeric(.x)))) %>% 
  mutate(path1 = paste0(sprintf('%0.3f', fit_g1), " (", sprintf('%0.3f', lcb_g1), ", ", sprintf('%0.3f', ucb_g1), ")"), 
         path2 = paste0(sprintf('%0.3f', fit_g2), " (", sprintf('%0.3f', lcb_g2), ", ", sprintf('%0.3f', ucb_g2), ")")) %>% 
  select(type, g1, g2, path1, path2) %>% 
  arrange(type, g1, g2)

write_csv(dia, file.path("/ihome/spurkayastha/soumik/2025_biometrika/output/element_2/bpd_clean.csv"))
write_csv(sys, file.path("/ihome/spurkayastha/soumik/2025_biometrika/output/element_2/bps_clean.csv"))
