rm(list = ls())
pacman::p_load(tidyverse, corrplot, reshape2, ggsci, patchwork, latex2exp, bnlearn)

source("/ihome/spurkayastha/soumik/2025_biometrika/code/helper.R")

bp <-  read_csv(file.path("/ihome/spurkayastha/soumik/2025_biometrika/data/outcome_bp.csv")) %>% rename(pid = foliocc) 
dia <- read_csv(file.path("/ihome/spurkayastha/soumik/2025_biometrika/data/element_dia.csv"))
sys <- read_csv(file.path("/ihome/spurkayastha/soumik/2025_biometrika/data/element_sys.csv"))

dia <- dia %>% pivot_wider(values_from = m, names_from = gene) %>% 
  inner_join(bp %>% dplyr::select(pid, bpd))

parms <- expand.grid(1:8, 1:8)
parms <- parms %>% filter(Var1 < Var2)

output <- tibble(x1 = "", x2 = "", x3 = "", x4 = "")

op_text_con <- NULL
op_text_exp <- NULL
op_text <- NULL


for(i in 1:nrow(parms)){
  gene1 = parms[i, 1]
  gene2 = parms[i, 2]
  
  x = dia %>% pull((1 + gene1))
  y = dia %>% pull((1 + gene2))
  bpd = dia %>% pull(bpd)
  
  x = affine(rank(x)/length(x))
  y = affine(rank(y)/length(y))
  bpd = affine(rank(bpd)/length(bpd))
  
  op_bpd <- el_overall2(cbind(x, y, bpd), neval = 50)
  
  op_text_con <- c(op_text_con, 
                   c(paste0("bpd_", colnames(dia)[1+gene1], "_", colnames(dia)[1+gene2]), 
                     score_helper(el_wrapper(op_bpd, "con"))))
  
  op_text_exp <- c(op_text_exp, 
                   c(paste0("bpd_", colnames(dia)[1+gene1], "_", colnames(dia)[1+gene2]), 
                     score_helper(el_wrapper(op_bpd, "exp"))))
  
  filename1 <- file.path("/ihome/spurkayastha/soumik/2025_biometrika/output/element_2/diastolic", 
                         paste0("bpd_con_", colnames(dia)[1+gene1], "_", colnames(dia)[1+gene2], ".csv"))
  
  filename2 <- file.path("/ihome/spurkayastha/soumik/2025_biometrika/output/element_2/diastolic", 
                         paste0("bpd_exp_", colnames(dia)[1+gene1], "_", colnames(dia)[1+gene2], ".csv"))
  
  filename3 <- file.path("/ihome/spurkayastha/soumik/2025_biometrika/output/element_2/diastolic", 
                         paste0("image_bpd_", colnames(dia)[1+gene1], "_", colnames(dia)[1+gene2], ".pdf"))
  
  op_bpd[[2]] <- op_bpd[[1]] %>% 
    rowwise() %>% 
    mutate(conditioner = ifelse(conditioner == "Y", 
                                paste0(colnames(dia)[1+gene1], " \u2192 ", "DBP | ", colnames(dia)[1+gene2]), 
                                paste0(colnames(dia)[1+gene2], " \u2192 ", "DBP | ", colnames(dia)[1+gene1]))) %>% 
    ggplot(aes(x = grid, group = conditioner)) + 
    #geom_point(aes(y = estim, color = conditioner), alpha = 0.5) + 
    geom_point(aes(y = fitted, color = conditioner)) + 
    geom_line(aes(y = fitted)) + 
    geom_ribbon(aes(ymin = lcb, ymax = ucb, fill = conditioner), alpha = 0.5) + 
    facet_grid(cols = vars(conditioner)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    scale_color_aaas() + 
    scale_fill_aaas() + 
    theme_bw() + 
    theme(legend.position = "none", 
          strip.background = element_rect(fill = "black"), 
          strip.text = element_text(color = "white", face = "bold", size = 12), 
          axis.title = element_text(face = "bold", size = 14),
          axis.text = element_text(face = "bold", size = 14)) + 
    labs(x = "", y = "Conditional asymmetry coefficient")
  
  suppressWarnings(ggsave(filename3, op_bpd[[2]]))
  
  temp_c <- op_bpd[[1]] %>% group_by(conditioner) %>% filter(estim == min(estim)) %>% ungroup() %>% 
    select(c(fitted, lcb, ucb))
  
  temp_e <- op_bpd[[1]] %>% group_by(conditioner) %>% filter(estim == max(estim)) %>% ungroup() %>% 
    select(c(fitted, lcb, ucb))
  
  op_text <- rbind(op_text, 
  rbind(as_tibble(matrix(unlist(c(colnames(dia)[1+gene1], colnames(dia)[1+gene2], cbind(temp_c[1,], temp_c[2,]), score_helper(el_wrapper(op_bpd, "con")))), byrow = F, ncol = 9)) %>% 
    add_column(type = "c"), 
    as_tibble(matrix(unlist(c(colnames(dia)[1+gene1], colnames(dia)[1+gene2], cbind(temp_e[1,], temp_e[2,]), score_helper(el_wrapper(op_bpd, "exp")))), byrow = F, ncol = 9)) %>% 
      add_column(type = "e")) %>% 
    rename(g1 = V1, g2 = V2, 
           fit_g1 = V3, lcb_g1 = V4, ucb_g1 = V5, 
           fit_g2 = V6, lcb_g2 = V7, ucb_g2 = V8, 
           ind = V9) %>% 
    select(g1, g2, type, ind, fit_g1, lcb_g1, ucb_g1, fit_g2, lcb_g2, ucb_g2))
}

write_csv(data.frame(matrix(op_text_con, byrow = T, ncol = 2)), file.path("/ihome/spurkayastha/soumik/2025_biometrika/output/element_2/bpd_con_summary.csv"))
write_csv(data.frame(matrix(op_text_exp, byrow = T, ncol = 2)), file.path("/ihome/spurkayastha/soumik/2025_biometrika/output/element_2/bpd_exp_summary.csv"))
write_csv(data.frame(op_text %>% filter(ind == 1) %>% mutate(across(contains("_"), ~as.numeric(.x)))), 
          file.path("/ihome/spurkayastha/soumik/2025_biometrika/output/element_2/bpd_summary.csv"))
