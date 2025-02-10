pacman::p_load(CAM, KernSmooth, here, rmi, tidyverse, bnlearn, pracma, igraph)
source("/ihome/spurkayastha/soumik/2025_biometrika/code/helper.R")

parms = expand.grid(n = c(100, 500, 1000), 
                    rho =c(0, 0.05, 0.10), 
                    i= seq(1:6), j=seq(1:6))


# i <- ifelse(is.na(as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))),
#             sample(1:nrow(parms), 1),
#             as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")))

i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

#i = 5
#i = 59
#i = 113

n = parms[i, 1]
rho = parms[i, 2]
pat_xz = parms[i, 3]
pat_yz = parms[i, 4]

op <- NULL

for(j in 1:250){
  input <- data_gen_step1(n, rho)
  output <- affine(data_gen_step2(input[,1], pat_xz) + data_gen_step2(input[,2], pat_yz))
  data <- cbind(input, output)
  colnames(data) <- c("x", "y", "z")
  
  op_el <- el_overall2(data, neval = 100) ## our method
  op_pc <- pc.stable(data.frame(data)) ### constraint-based
  op_hc <- hc(data.frame(data)) ### score-based
  op_rs <- rsmax2(data.frame(data)) ## hybrid
  
  op <- rbind(op, 
              c(graph_sim_calc(graph_adj_maker(el_wrapper(op_el))), 
                graph_sim_calc(graph_adj_maker(op_pc$arcs)), 
                graph_sim_calc(graph_adj_maker(op_hc$arcs)), 
                graph_sim_calc(graph_adj_maker(op_rs$arcs))))
}

metrics <- apply(data.frame((op)), 2, mean)
parmvals <- c(n, rho, pat_xz, pat_yz)

temp <- c(parmvals, metrics)

names(temp) <- c("n", "rho", "pat_xz", "pat_yz", 
               paste0(rep(c("EL", "PC", "HC", "RS"), each = 3), 
                         rep("_", times = 12), 
                         rep(c("Frob", "Hamm", "Edit"), times = 4)))

filename <- file.path("/ihome/spurkayastha/soumik/2025_biometrika/output/simulation", 
                     paste0("sim_", i, ".csv"))

#t(data.frame(temp))[c(1,2,3,4,5,8,11,14)]
write_csv((data.frame(temp)), filename)


