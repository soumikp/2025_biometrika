affine <- function(x){
  return((x - min(x))/(max(x) - min(x)))
}

entropy_estim <- function(x){
  return(mean(-log(x)*x))
}

entcalc <- function(p){
  p <- p[p > 0]
  return(mean(-log(p)*p))
}

el_estim <- function(data, frac = 10){
  
  x <- data[,1]
  y <- data[,2]
  z <- data[,3]
  
  grid <- seq(0, 1, length.out = round(length(x)/frac))
  
  c_xz_y <- NULL
  c_yz_x <- NULL
  #c_xy_z <- NULL
  
  for(i in 1:length(grid)){
    y_indices <- which(abs(y - grid[i]) <= 0.10)
    c_xz_y <- c(c_xz_y, rmi::lnn_entropy(x[y_indices]) - rmi::lnn_entropy(z[y_indices]))
    
    x_indices <- which(abs(x - grid[i]) <= 0.10)
    c_yz_x <- c(c_yz_x, rmi::lnn_entropy(y[x_indices]) - rmi::lnn_entropy(z[x_indices]))
    
    #z_indices <- which(abs(z - grid[i]) <= 0.25)
    #c_xy_z <- c(c_xy_z, rmi::lnn_entropy(x[z_indices]) - rmi::lnn_entropy(y[z_indices]))
  }
  return(cbind(grid, c_xz_y, c_yz_x#, c_xy_z
  ))
}

fit_lq <- function(x, y, h){
  l <- length(x)
  beta0 <- rep(0, l)
  beta1 <- rep(0, l)
  for(i in 1:l) {
    x.reg <- x - x[i]
    w <- dnorm(x - x[i], 0, h)
    fit <- lm(y ~ x.reg + I(x.reg^ 2), weights = w)
    beta0[i] <- fit$coe[1]
    beta1[i] <- fit$coe[2]
  }
  beta <- cbind(beta0, beta1)
  return(beta)
}

add_ci <- function(x, y, h, beta){
  l <- length(x)
  beta0 <- beta[, 1]
  beta1 <- beta[, 2]
  ##Test for equilibruim points
  w <- rep( 0, l)
  diag <- rep( 0, l)
  
  upperCI_0 <- rep( 0, l)
  lowerCI_0 <- rep( 0, l)
  
  upperCI_1 <- rep( 0, l)
  lowerCI_1 <- rep( 0, l)
  
  se_0 <- rep( 0, l)
  se_1 <- rep( 0, l)
  
  z_0 <- rep( 0, l)
  z_1 <- rep( 0, l)
  
  p_0 <- rep( 0 ,l)
  p_1 <- rep( 0 ,l)
  options( object.size = 1000000000)
  ##Estimate sigma^2
  for( i in 1:l){
    Xi <- cbind( rep( 1, l), x - x[i], ( x - x[i]) ^ 2)
    ker <- dnorm( Xi[,2], 0, h)
    A <- matrix( 0, ncol = 3, nrow = 3)
    A[1,1] <- sum( ker)
    A[1,2] <- ker %*% Xi[,2]
    A[2,1] <- A[1,2]
    A[1,3] <- ker %*% Xi[,3]
    A[2,2] <- A[1,3]
    A[3,1] <- A[1,3]
    A[2,3] <- ker %*% Xi[,2] ^ 3
    A[3,2] <- A[2,3]
    A[3,3] <- ker %*% Xi[,3] ^ 2
    B <- solve( A)[1,]
    C <- rbind( ker, ker * Xi[,2], ker * Xi[,3])
    wi <- B %*% C
    diag[i] <- wi[i]
    w <- rbind( w, wi)
  }
  w <- w[ 2:( l + 1), ]
  second <- sum( w ^ 2)
  first <- 2 * sum( diag)
  v <- first - second
  vari <- 1 / ( l - v ) * sum( ( y - beta0) ^ 2)
  ##Calculate the 95% confidence band
  for( i in 1:l) {
    X <- cbind( rep( 1, l), x - x[i], ( x - x[i]) ^ 2)
    kernel <- dnorm( X[, 2], 0, h)
    An <- matrix( 0, ncol = 3, nrow = 3)
    Bn <- matrix( 0, ncol = 3, nrow = 3)
    An[1,1] <- sum( kernel) / l
    An[1,2] <- kernel %*% X[,2] / l
    An[2,1] <- An[1,2]
    An[1,3] <- kernel %*% X[,3] / l
    An[2,2] <- An[1,3]
    An[3,1] <- An[1,3]
    An[2,3] <- kernel %*% X[,2] ^ 3 / l
    An[3,2] <- An[2,3]
    An[3,3] <- kernel %*% X[,3] ^ 2 / l
    kernel2 <- kernel ^ 2
    Bn[1,1] <- sum( kernel2) / l / l
    Bn[1,2] <- kernel2 %*% X[,2] / l / l
    Bn[2,1] <- Bn[1,2]
    Bn[1,3] <- kernel2 %*% X[,3] / l / l
    Bn[2,2] <- Bn[1,3]
    Bn[3,1] <- Bn[1,3]
    Bn[2,3] <- kernel2 %*% X[,2] ^ 3 / l / l
    Bn[3,2] <- Bn[2,3]
    Bn[3,3] <- kernel2 %*% X[,3] ^ 2 / l / l
    sol <- solve( An)
    temp <- sol %*% Bn %*% sol
    
    temp1 <- temp[1,1]
    se_0[i] <- sqrt( vari * temp1)
    z_0[i] <- beta0[i] / se_0[i]
    p_0[i] <- (1 - pnorm(z_0[i]))
    upperCI_0[i] <- beta0[i] + 1.96 * se_0[i]
    lowerCI_0[i] <- beta0[i] - 1.96 * se_0[i]
    
    temp2 <- temp[2,2]
    se_1[i] <- sqrt( vari * temp2)
    z_1[i] <- abs( beta1[i] / se_1[i])
    p_1[i] <- ( 1 - pnorm( z_1[i])) * 2
    upperCI_1[i] <- beta1[i] + qnorm(1-0.025/2) * se_1[i]
    lowerCI_1[i] <- beta1[i] - qnorm(1-0.025/2) * se_1[i]
    
    
  }
  upperCI_0 <- round(upperCI_0, 5)
  lowerCI_0 <- round(lowerCI_0, 5)
  
  upperCI_1 <- round(upperCI_1, 5)
  lowerCI_1 <- round(lowerCI_1, 5)
  
  p_0 <- round(p_0, 5)
  p_1 <- round(p_1, 5)
  
  CIp_0 <- cbind( upperCI_0, lowerCI_0, p_0)
  CIp_1 <- cbind( upperCI_1, lowerCI_1, p_1)
  return(cbind(CIp_0, CIp_1))
}

find_bw <- function(x, y){
  return(nprobust::lpbwselect(x = x, y = y,
                              eval = x[y == min(y)])$bws[2])
}

data_gen_step1 <- function(n, rho){
  data <- MASS::mvrnorm(n, c(0, 0), cbind(c(1, rho), c(rho, 1)))
  x <- rank(data[,1])/(n+1)
  y <- rank(data[,2])/(n+1)
  return(cbind(x, y))
}

data_gen_step2 <- function(x, type = 1){
  if(type == 1){
    y <- sqrt(x)
  }else if(type == 2){
    y <- x^2
  }else if(type == 3){
    y <- (2*x - 1)/(x^2 - x + 5)
  }else if(type == 4){
    y <- log(1/(1 + exp(-x)))
  }else if(type == 5){
    y <- (sin(pi*(x - 0.5)) + exp(x))/(1 + exp(x))
  }else if(type == 6){
    y <- (x^2 + sqrt(x))/(1 + sqrt(1 - x^2))
  }else if(type == 7){
    y <- exp(x)
  }else if(type == 8){
    y <- sin(pi*(x - 0.5))
  }
  return(y)
}

el_overall <- function(data, frac = 10){
  
  data2 <- data
  
  step1 <- el_estim(data2, frac = frac)
  
  step2a <- fit_lq(step1[,1], step1[,2],  find_bw(step1[,1], step1[,2]))
  step2b <- fit_lq(step1[,1], step1[,3],  find_bw(step1[,1], step1[,3]))
  #step2c <- fit_lq(step1[,1], step1[,4],  find_bw(step1[,1], step1[,4]))
  
  step3a <- add_ci(step1[,1], step1[,2],  find_bw(step1[,1], step1[,2]), step2a)
  step3b <- add_ci(step1[,1], step1[,3],  find_bw(step1[,1], step1[,3]), step2b)
  #step3c <- add_ci(step1[,1], step1[,4],  find_bw(step1[,1], step1[,4]), step2c)
  
  
  op_el <- as_tibble(rbind(cbind(step1[,c(1, 2)], step2a[,1], step3a[,c(1, 2)]),
                           cbind(step1[,c(1, 3)], step2b[,1], step3b[,c(1, 2)])#,
                           #cbind(step1[,c(1, 4)], step2c[,1], step3c[,c(1, 2)])
  )) %>% 
    rename(estim = c_xz_y, 
           fitted = V3, 
           lcb = lowerCI_0, 
           ucb = upperCI_0) %>% 
    add_column(conditioner = rep(c("Y", "X"#, "Z"
    ), each = nrow(step1)))
  
  plot_el <- op_el %>% 
    ggplot(aes(x = grid, group = conditioner)) + 
    geom_point(aes(y = estim, color = conditioner), alpha = 0.5) + 
    geom_point(aes(y = fitted, color = conditioner)) + 
    geom_line(aes(y = fitted)) + 
    geom_ribbon(aes(ymin = lcb, ymax = ucb, fill = conditioner), alpha = 0.7) + 
    facet_grid(cols = vars(conditioner)) + 
    geom_hline(yintercept = 0, linetype = "dashed")
  
  return(list(op_el, plot_el))
}

el_overall2 <- function(data, neval = 100){
  
  cond_y <- apply(hdrcde::cde(data[,2], data[,1], nxmargin = neval, nymargin = 100)$z, 1, entcalc) - apply(hdrcde::cde(data[,2], data[,3], nxmargin = neval, nymargin = 100)$z, 1, entcalc)
  grid_y <- hdrcde::cde(data[,2], data[,1], nxmargin = neval, nymargin = 100)$x
  cond_x <- apply(hdrcde::cde(data[,1], data[,2], nxmargin = neval, nymargin = 100)$z, 1, entcalc) - apply(hdrcde::cde(data[,1], data[,3], nxmargin = neval, nymargin = 100)$z, 1, entcalc)
  grid_x <- hdrcde::cde(data[,1], data[,2], nxmargin = neval, nymargin = 100)$x
  
  step2a <- fit_lq(grid_y, cond_y,  find_bw(grid_y, cond_y))
  step2b <- fit_lq(grid_x, cond_x,  find_bw(grid_x, cond_x))
  #step2c <- fit_lq(step1[,1], step1[,4],  find_bw(step1[,1], step1[,4]))
  
  step3a <- add_ci(grid_y, cond_y,  find_bw(grid_y, cond_y), step2a)
  step3b <- add_ci(grid_x, cond_x,  find_bw(grid_x, cond_x), step2b)
  #step3c <- add_ci(step1[,1], step1[,4],  find_bw(step1[,1], step1[,4]), step2c)
  
  
  op_el <- as_tibble(rbind(cbind(grid_y, cond_y, step2a[,1], step3a[,c(1, 2)]),
                           cbind(grid_x, cond_x, step2b[,1], step3b[,c(1, 2)])#,
                           #cbind(step1[,c(1, 4)], step2c[,1], step3c[,c(1, 2)])
  )) %>% 
    rename(estim = cond_y, 
           grid = grid_y, 
           fitted = V3, 
           lcb = lowerCI_0, 
           ucb = upperCI_0) %>% 
    add_column(conditioner = rep(c("Y", "X"#, "Z"
    ), each = nrow(step2a)))
  
  plot_el <- op_el %>% 
    ggplot(aes(x = grid, group = conditioner)) + 
    geom_point(aes(y = estim, color = conditioner), alpha = 0.5) + 
    geom_point(aes(y = fitted, color = conditioner)) + 
    geom_line(aes(y = fitted)) + 
    geom_ribbon(aes(ymin = lcb, ymax = ucb, fill = conditioner), alpha = 0.7) + 
    facet_grid(cols = vars(conditioner), scales = "free_x") + 
    geom_hline(yintercept = 0, linetype = "dashed")
  
  return(list(op_el, plot_el))
}

el_wrapper <- function(op, funtype = "con"){
  
  if(funtype == "con"){
    temp <- op[[1]] %>% mutate(sign = ifelse(lcb > 0, 1, 0)) %>% group_by(conditioner) %>% filter(estim == min(estim)) %>% filter(sign == 1)
  }else{
    temp <- op[[1]] %>% mutate(sign = ifelse(ucb < 0, 1, 0)) %>% group_by(conditioner) %>% filter(estim == max(estim)) %>% filter(sign == 1)
  }
  
  op <- NULL
  
  if(nrow(temp) != 0){
    if("Y" %in% temp$conditioner){
      op <- rbind(op, c("x", "z"))
    }
    
    if("X" %in% temp$conditioner){
      op <- rbind(op, c("y", "z"))
    }
    
    if("Z" %in% temp$conditioner){
      op <- rbind(op, c("x", "y"))
    }
    
    colnames(op) <- c("from", "to")
    
  }
  
  
  return(op)
}

score_helper <- function(data){
  if(is.null(data))
    return(0)
  if(nrow(data) != 2){
    return(0)
  }else if(sum(sort(data[,1]) == c("x", "y")) != 2){
    return(0)
  }else if(sum(sort(data[,2]) == c("z", "z")) != 2){
    return(0)
  }else{
    return(1)
  }
}

graph_adj_maker <- function(op){
  to <- from <- c("x", "y", "z")
  adj <- matrix(0, ncol = 3, nrow = 3)
  adj[cbind(match(op[,1], from), match(op[,2], to))] <- 1
  return(adj)
}

graph_sim_calc <- function(mat){
  ## mat is the input adjacency matrix. 
  
  target <- matrix(0, ncol = 3, nrow = 3)
  target[1,3] <- target[2,3] <- 1
  
  
  frobenius <- sqrt(sum(mat-target)^2)
  hamming <- sum(mat!=target)
  edit <- sum(abs(mat-target))
  
  return(c(frobenius, hamming, edit))
}


