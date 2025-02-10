pacman::p_load(purrr, tidyverse, data.table, here)

files <- list.files(file.path(here(), "2025_biometrika/output/simulation"), 
                    full.names = TRUE)

op <- lapply(files, purrr::map_dfr, fread)

helper <- function(data){
  return(unlist(c(data[1, 1:4],
    data %>% summarise(el = mean(X6), 
                     pc = mean(X7), 
                     hc = mean(X8), 
                     rsm = mean(X9)))))
}


op2 <- as_tibble(matrix(unlist(op), byrow = T, ncol= 16)) %>% 
  select(c(1, 2, 3, 4, 5, 8, 11))

helper_text <- function(x){
  return(paste0(x[1], " (", x[2], ") [", x[3], "]"))
}

op3 <- op2 %>% 
  rename(n = V1, r = V2, xz = V3, yz = V4, el = V5, pc = V8, hc = V11) %>% 
  arrange(xz, yz) %>% 
  mutate(el_text = sprintf("%.2f", el),
         pc_text = sprintf("%.2f", pc),
         hc_text = sprintf("%.2f", hc)) %>% 
  select(-c(el, pc, hc)) %>% 
  group_by(n, xz, yz) %>% 
  summarise(el = str_c(str_c(el_text[1], " (", el_text[2], ") [", el_text[3], "]"), collapse = ""),
            hc = str_c(str_c(hc_text[1], " (", hc_text[2], ") [", hc_text[3], "]"), collapse = ""),
            pc = str_c(str_c(pc_text[1], " (", pc_text[2], ") [", pc_text[3], "]"), collapse = "")) %>% 
  ungroup() %>% 
  arrange(xz, yz, n)

write_csv(op3 %>% filter(n == 100), 
          file.path(here(), "2025_biometrika/output/simulation_processed_100.csv"))

write_csv(op3 %>% filter(n == 500), 
          file.path(here(), "2025_biometrika/output/simulation_processed_500.csv"))

write_csv(op3 %>% filter(n == 1000), 
          file.path(here(), "2025_biometrika/output/simulation_processed_1000.csv"))



  
 