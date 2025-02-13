#===============================================================
# Author: Andreas Anastasiou
# Project: Non-parametric multiple change-point detection

# Description: The code below provides the NPID routine as well as some main
# functions that are necessary in order to be able to run the simulation study
# in the paper.

#===============================================================
rm(list = ls())   # Clear environment
library(stats)
library(IDetect)

#==============================
# Some utility functions
#==============================

## A function used in case we want to rescale the obtained values based on the
## standard deviation
sigma_p <- function(x){
  p <- mean(x)
  res <- ifelse(((p >= 0.1) & (p <= 0.9)), sqrt(p * (1-p)), 0.3)
  res
}

## Creating the binary sequences
Binary_seq  <- function(x, grid.size = NULL){
  lx <- length(x)
  mx  <- min(x)
  Mx <- max(x)
  sx <- sort(x, index = TRUE)
  cts <- cumsum(rle(sx[[1]])[[1]])
  lsx <- length(unique(sx[[1]]))
  if (lsx == 1){
    stop("The given data sequence should not have all values the same")}
  else{
    if (is.null(grid.size)){
      ber_grid_na <- 0.5*(sx[[1]][1:lsx] + sx[[1]][2:(lsx+1)])
      ber_grid <- ber_grid_na[-lsx]
      result <- matrix(1, lx, lsx - 1)
      for (i in 1:(lsx-1)){
        result[sx[[2]][1:cts[i]],i] <- 0
      }
    }
    else{ber_grid <- seq(mx + (Mx - mx)/(grid.size+1), Mx - (Mx - mx)/(grid.size+1),(Mx - mx)/(grid.size+1))
    result <- matrix(1, lx, grid.size)
    i <- 1
    while(i <= grid.size){
      result[which(x < ber_grid[i]),i] <- 0
      i <- i+1
    }
    }
    return(list(result, ber_grid))
  }
}


## Calculating cumulative sums for the elements of a given matrix
cumul_sum_matrix <- function(x){
  lx <- dim(x)[1]
  l_grid <- dim(x)[2]
  cumul_sum <- matrix(NA, lx, l_grid)
  i <- 1
  while(i <= l_grid){
    cumul_sum[,i] <- cumsum(x[,i])
    i <- i+1
  }
  return(cumul_sum)
}

## Calculating the cusum function values regarding all the columns of the input xber,
## which consist of binary sequences.
cusum_matrix <- function(xber, y = cumul_sum_matrix(xber), sp=1, ep = dim(y)[1], resc = FALSE){
  l_seq <- ep - sp + 1
  l_grid <- dim(y)[2]
  r1 <- ((l_seq - 1) : 1)/ rev(((l_seq - 1) : 1))
  q1 <- r1 / l_seq
  q2 = r1 * l_seq
  b <- sp:(ep - 1)
  result <- matrix(NA, l_seq - 1, l_grid)
  if (sp == 1){
    if (resc == FALSE) {
      for (i in 1:l_grid){
        result[,i] <- abs(sqrt( q1) * y[b, i] - sqrt( 1 / q2) * (y[ep, i] - y[b, i]))}
    }
    else{
      for (i in 1:l_grid){
        result[,i] <- abs(sqrt( q1) * (y[b, i]) - sqrt( 1 / q2) * (y[ep, i] - y[b, i]))/ sigma_p(xber[1:ep,i])}
    }
  }
  else{
    if (resc == FALSE) {
      for (i in 1:l_grid){
        result[,i] <- abs(sqrt( q1) * (y[b, i] - y[sp-1, i]) - sqrt( 1 / q2) * (y[ep, i] - y[b, i]))}
    }
    else{
      for (i in 1:l_grid){
        result[,i] <- abs(sqrt( q1) * (y[b, i] - y[sp-1, i]) - sqrt( 1 / q2) * (y[ep, i] - y[b, i]))/ sigma_p(xber[sp:ep,i])}
    }
  }
  return(result)
}

## This function finds the value of the relevant non-parametric CUSUM at b, when
## working in the interval [s,e)
CUSUM_at_point <- function(x, y = cumul_sum_matrix(x), s, e, b, rescale = FALSE) {
  if (!(is.numeric(x))){
    stop("The input in `x' should be a numeric matrix with values either equal to 0 or equal to 1.")
  }
  if ( (length(s) != length(b)) || (length(s) != length(e)) || (length(e) != length(b))){
    stop("The vectors s, b, e, should be of the same length")
  }
  if (any(s < 1) | any(b < 1) | any(e < 1)){
    stop("The entries of the vectors s, b, e should be positive integers.")
  }
  if (any(s > b) | any(b >= e)){
    stop("The value for b should be in the interval [s,e)")
  }
  if ( (any(abs( (s - round(s))) > .Machine$double.eps ^ 0.5))
       || (any(abs( (b - round(b))) > .Machine$double.eps ^ 0.5))
       || (any(abs( (e - round(e))) > .Machine$double.eps ^ 0.5))){
    stop("The input values  for s, b, and  e should be positive integers.")
  }
  l_grid <- dim(y)[2]
  l <- numeric()
  alpha <- numeric()
  beta <- numeric()
  result_matrix <- matrix(NA, length(b), l_grid)
  d <- matrix(NA, length(b), l_grid)
  for (j in 1:length(b)) {
    l[j] <- e[j] - s[j] + 1
    alpha[j] <- b[j] - s[j] + 1
    beta[j] <- e[j] - b[j]
    if(s[j] == 1){
      d[j,] <- rep(0, l_grid)}
    else{
      d[j,] <- y[s[j]-1, ]}
    result_matrix[j,1:l_grid] <- abs(sqrt( beta[j] / (l[j] * alpha[j])) * (y[b[j], 1:l_grid] - d[j, 1:l_grid]) - sqrt( alpha[j] / (l[j] * beta[j])) * (y[e[j], 1:l_grid] - y[b[j], 1:l_grid]))
    if (rescale == TRUE) {
      for (i in 1:l_grid){result_matrix[j,i] <- result_matrix[j,i]/ sigma_p(x[s[j]:e[j],i])}
    } 
  }
  return(apply(result_matrix, 1, max))
}

## The maximum joint log-likelihood for a given set of change-points on given data.
loglik_values <- function(x, b1 = Binary_seq(x)[[1]], y = cumul_sum_matrix(b1), cpt,
                          no.cpt = length(cpt)){
  lx<- length(x)
  if (no.cpt == 0){
    hat_Fkl <- (1 - (y[lx,2:(dim(y)[2])])/(lx))
    result_matr <- (lx-1) * (xlogx(hat_Fkl) + xlogx(1 - hat_Fkl))/((1:(ncol(y) - 1))*((ncol(y)):2))}
  else{
    crucial_points <- c(1,sort(cpt),lx)
    dif_cruc <- diff(crucial_points)
    hat_Fkl <- matrix(NA, no.cpt + 1, ncol(y) - 1)
    result_matr <- matrix(NA, no.cpt + 1, ncol(y) - 1)
    for (l in 1:(ncol(y) - 1)){
      hat_Fkl[,l] <- (1 - (y[crucial_points[2:(no.cpt + 2)],l+1] - c(0,y[crucial_points[2:(no.cpt+1)] - 1,l+1]))/(dif_cruc + 1))
      result_matr[,l] <- dif_cruc[1:(no.cpt+1)] * (xlogx(hat_Fkl[,l]) + xlogx(1 - hat_Fkl[,l]))/(l*(lx - l))}
  }
  result <- lx * sum(sum(result_matr))
  result
}

#=========================================================
# Main functions for the thresholding-based decision rule
#=========================================================

## The approach based on the L_inf mean-dominant norm
L_inf_norm <- function(x, thr_const = 0.9, thr_fin = thr_const * sqrt(log(length(x))),
                           s = 1, e = length(x), points = 10, r_e_points = seq(points, l, points),
                           l_e_points = seq(l - points + 1, 1, -points), k_l = 1, k_r = 1, gs= NULL,
                           rescale = FALSE, brn_list = Binary_seq(x, grid.size = gs)[[1]], 
                           cs = cumul_sum_matrix(brn_list)) {
  if (!(is.numeric(x))){
    stop("The input in `x' should be a numeric vector containing the data in
         which you would like to find change-points.")
  }
  if ( (thr_const <= 0) || (points <= 0)){
    stop("The threshold constant as well as the `points' argument that represents the
         magnitude of the expansion for the intervals should be positive numbers.")
  }
  if (abs(points - round(points)) > .Machine$double.eps ^ 0.5){
    warning("The input for `points' should be a positive integer. If it is a positive real
            number then the integer part of the given number is used as the value of `points'.")
  }
  if (e - s <= 1) {
    cpt <- 0
  }
  else {
    grid.size <- ncol(brn_list)
    l <- length(x)
    chp <- 0
    moving_points <- IDetect::s_e_points(r_e_points, l_e_points, s, e)
    right_points <- moving_points[[1]]
    left_points <- moving_points[[2]]
    lur <- length(left_points)
    rur <- length(right_points)
    if (k_r < k_l) {
      while ( (chp == 0) & (k_r < min(k_l, rur))) {
        ipcr <- cusum_matrix(x = brn_list, y = cs, sp = s, ep = right_points[k_r], resc = rescale)
        sipcr <- apply(ipcr, 1, max)
        pos_r <- which.max(sipcr) + s - 1
        CUSUM_r <- sipcr[pos_r - s + 1]
        if (CUSUM_r > thr_fin) {
          chp <- pos_r
        } else {
          k_r <- k_r + 1
        }
      }
    }
    if (k_l < k_r) {
      while ( (chp == 0) & (k_l < min(k_r, lur))) {
        ipcl <- cusum_matrix(x = brn_list, y = cs, sp = left_points[k_l], ep = e, resc = rescale)
        sipcl <- apply(ipcl, 1, max)
        pos_l <- which.max(sipcl) + left_points[k_l] - 1
        CUSUM_l <- sipcl[pos_l - left_points[k_l] + 1]
        if (CUSUM_l > thr_fin) {
          chp <- pos_l
        } else {
          k_l <- k_l + 1
        }
      }
    }
    if (chp == 0) {
      while ( (chp == 0) & (k_l <= lur) & (k_r <= rur)) {
        ipcr <- cusum_matrix(x = brn_list, y = cs, sp = s, ep = right_points[k_r], resc = rescale)
        sipcr <- apply(ipcr, 1, max)
        pos_r <- which.max(sipcr) + s - 1
        CUSUM_r <- sipcr[pos_r - s + 1]
        if (CUSUM_r > thr_fin) {
          chp <- pos_r
        } else {
          ipcl <- cusum_matrix(x = brn_list, y = cs, sp = left_points[k_l], ep = e, resc = rescale)
          sipcl <- apply(ipcl, 1, max)
          pos_l <- which.max(sipcl) + left_points[k_l] - 1
          CUSUM_l <- sipcl[pos_l - left_points[k_l] + 1]
          if (CUSUM_l > thr_fin) {
            chp <- pos_l
          } else {
            k_r <- k_r + 1
            k_l <- k_l + 1
          }
        }
      }
    }
    if (chp != 0) {
      if (chp > ( (e + s) / 2)) {
        r <- L_inf_norm(x, s = s, e = chp, points = points, r_e_points = r_e_points, l_e_points = l_e_points,
                            thr_fin = thr_fin, k_r = k_r, k_l = 1,brn_list = brn_list, cs = cs,
                            rescale = rescale)[[1]]
      } else {
        r <- L_inf_norm(x, s = chp + 1, e = e, points = points, r_e_points = r_e_points, l_e_points = l_e_points,
                            thr_fin = thr_fin, k_r = 1, k_l = max(1, k_l - 1), brn_list = brn_list, cs = cs, 
                            rescale = rescale)[[1]]
      }
      cpt <- c(chp, r)
    } else {
      cpt <- chp
    }
  }
  cpt <- cpt[cpt != 0]
  return(list(changepoints= sort(cpt), bernoulli_sequence = brn_list, cms_matrix = cs))
}

## The approach based on the L_2 mean-dominant norm
L_2_norm <- function(x, thr_const = 0.6, thr_fin = thr_const * sqrt(log(length(x))),
                            s = 1, e = length(x), points = 10, r_e_points = seq(points, l, points),
                            l_e_points = seq(l - points + 1, 1, -points), k_l = 1, k_r = 1, gs= NULL,
                            rescale = FALSE, brn_list = Binary_seq(x, grid.size = gs)[[1]], 
                            cs = cumul_sum_matrix(brn_list)) {
  if (!(is.numeric(x))){
    stop("The input in `x' should be a numeric vector containing the data in
         which you would like to find change-points.")
  }
  if ( (thr_const <= 0) || (points <= 0)){
    stop("The threshold constant as well as the `points' argument that represents the
         magnitude of the expansion for the intervals should be positive numbers.")
  }
  if (abs(points - round(points)) > .Machine$double.eps ^ 0.5){
    warning("The input for `points' should be a positive integer. If it is a positive real
            number then the integer part of the given number is used as the value of `points'.")
  }
  if (e - s <= 1) {
    cpt <- 0
  }
  else {
    grid.size <- ncol(brn_list)
    l <- length(x)
    chp <- 0
    moving_points <- IDetect::s_e_points(r_e_points, l_e_points, s, e)
    right_points <- moving_points[[1]]
    left_points <- moving_points[[2]]
    lur <- length(left_points)
    rur <- length(right_points)
    if (k_r < k_l) {
      while ( (chp == 0) & (k_r < min(k_l, rur))) {
        ipcr <- cusum_matrix(x = brn_list, y = cs, sp = s, ep = right_points[k_r], resc = rescale)
        sipcr <- sqrt(rowSums(ipcr^2))/sqrt(grid.size)
        pos_r <- which.max(sipcr) + s - 1
        CUSUM_r <- sipcr[pos_r - s + 1]
        if (CUSUM_r > thr_fin) {
          chp <- pos_r
        } else {
          k_r <- k_r + 1
        }
      }
    }
    if (k_l < k_r) {
      while ( (chp == 0) & (k_l < min(k_r, lur))) {
        ipcl <- cusum_matrix(x = brn_list, y = cs, sp = left_points[k_l], ep = e, resc = rescale)
        sipcl <- sqrt(rowSums(ipcl^2))/sqrt(grid.size)
        pos_l <- which.max(sipcl) + left_points[k_l] - 1
        CUSUM_l <- sipcl[pos_l - left_points[k_l] + 1]
        if (CUSUM_l > thr_fin) {
          chp <- pos_l
        } else {
          k_l <- k_l + 1
        }
      }
    }
    if (chp == 0) {
      while ( (chp == 0) & (k_l <= lur) & (k_r <= rur)) {
        ipcr <- cusum_matrix(x = brn_list, y = cs, sp = s, ep = right_points[k_r], resc = rescale)
        sipcr <- sqrt(rowSums(ipcr^2))/sqrt(grid.size)
        pos_r <- which.max(sipcr) + s - 1
        CUSUM_r <- sipcr[pos_r - s + 1]
        if (CUSUM_r > thr_fin) {
          chp <- pos_r
        } else {
          ipcl <- cusum_matrix(x = brn_list, y = cs, sp = left_points[k_l], ep = e, resc = rescale)
          sipcl <- sqrt(rowSums(ipcl^2))/sqrt(grid.size)
          pos_l <- which.max(sipcl) + left_points[k_l] - 1
          CUSUM_l <- sipcl[pos_l - left_points[k_l] + 1]
          if (CUSUM_l > thr_fin) {
            chp <- pos_l
          } else {
            k_r <- k_r + 1
            k_l <- k_l + 1
          }
        }
      }
    }
    if (chp != 0) {
      if (chp > ( (e + s) / 2)) {
        r <- L_2_norm(x, s = s, e = chp, points = points, r_e_points = r_e_points, l_e_points = l_e_points,
                             thr_fin = thr_fin, k_r = k_r, k_l = 1,brn_list = brn_list, cs = cs,
                             rescale = rescale)[[1]]
      } else {
        r <- L_2_norm(x, s = chp + 1, e = e, points = points, r_e_points = r_e_points, l_e_points = l_e_points,
                             thr_fin = thr_fin, k_r = 1, k_l = max(1, k_l - 1), brn_list = brn_list, cs = cs, 
                             rescale = rescale)[[1]]
      }
      cpt <- c(chp, r)
    } else {
      cpt <- chp
    }
  }
  cpt <- cpt[cpt != 0]
  return(list(changepoints= sort(cpt), bernoulli_sequence = brn_list, cms_matrix = cs))
}

## The L_inf window related variant in the case of very long data sequences
window_L_inf_norm <- function(xd, thr_con = NULL, c_win = 1000, w_points = 10, l_win = 5000, rsc = FALSE) {
  if (length(thr_con) == 0){thr_con = 0.8
  if (rsc == TRUE){thr_con = 1.9}}
  if (!(is.numeric(xd))){
    stop("The input in `xd' should be a numeric vector containing the data in
         which you would like to find change-points.")
  }
  if ( (thr_con <= 0) || (w_points <= 0) || (c_win <= 0) || (l_win <= 0)){
    stop("The threshold constant as well as the `w_points', `c_win', `l_win' arguments should
         be positive numbers.")
  }
  if ( (abs(w_points - round(w_points)) > .Machine$double.eps ^ 0.5)
       || (abs(c_win - round(c_win)) > .Machine$double.eps ^ 0.5)
       || (abs(l_win - round(l_win)) > .Machine$double.eps ^ 0.5)){
    warning("The input values  for `w_points', `c_win', and  `l_win' should be positive integers.
            If either of them is a positive real number then the integer part of the given number
            is used to obtain the result.")
  }
  lg <- length(xd)
  w_points <- as.integer(w_points)
  c_win <- min(lg, c_win)
  c_win <- as.integer(c_win)
  l_win <- as.integer(l_win)
  if (lg <= c_win) {
    u <- L_inf_norm(x = xd, thr_const = thr_con, points = w_points, rescale = rsc, gs = NULL)
    return(u)
  } else if ((lg > c_win) & (lg <= l_win)) {
    u <- L_inf_norm(x = xd, thr_const = thr_con, points = w_points, rescale = rsc, gs = c_win)
    return(u)
  }
  else {
    bs <- Binary_seq(xd, grid.size = c_win)[[1]]
    cs <- cumul_sum_matrix(bs)
    t <- thr_con * sqrt(log(lg))
    K <- ceiling(lg / c_win)
    tsm <- list()
    u <- list()
    ufin <- numeric()
    uaddition <- numeric(0)
    tsm[[1]] <- xd[1:c_win]
    ufin <- L_inf_norm(tsm[[1]], thr_fin = t, points = w_points, rescale = rsc)[[1]]
    if (CUSUM_at_point(x = bs,y = cs, s = max(1, c_win - (10 * w_points) + 1), e = min( (c_win + (10 * w_points)), lg), b = c_win) > t){
      uaddition <- c_win 
    }
    ufin <- c(ufin, uaddition)
    i <- 2
    while (i < K) {
      tsm[[i]] <- xd[( (i - 1) * c_win + 1):(i * c_win)]
      u[[i]] <- L_inf_norm(x = tsm[[i]], thr_fin = t, points = w_points, rescale = rsc)[[1]] + (i - 1) * c_win
      if (CUSUM_at_point(x = bs,y = cs, s = max(1, i * c_win - (10 * w_points) + 1), e = min(i * c_win + (10 * w_points), lg), b = i * c_win) > t){
        uaddition <- i*c_win
      }
      ufin <- c(ufin, u[[i]], uaddition)
      i <- i + 1
    }
    tsm[[K]] <- xd[( (K - 1) * c_win + 1):lg]
    u[[K]] <- L_inf_norm(tsm[[K]], thr_fin = t, points = w_points, rescale = rsc)[[1]] + (K - 1) * c_win
    ufinl <- c(ufin, u[[K]])
    return(list(cpt = sort(unique(ufinl)), bernoulli_sequence  = bs, cumul_sum = cs))
  }
}

## The L_2 window related variant in the case of very long data sequences
window_L_2_norm <- function(xd, thr_con = NULL, c_win = 1000, w_points = 10, l_win = 5000, rsc = FALSE) {
  if (length(thr_con) == 0){thr_con = 0.6
  if (rsc == TRUE){thr_con = 1}}
  if (!(is.numeric(xd))){
    stop("The input in `xd' should be a numeric vector containing the data in
         which you would like to find change-points.")
  }
  if ( (thr_con <= 0) || (w_points <= 0) || (c_win <= 0) || (l_win <= 0)){
    stop("The threshold constant as well as the `w_points', `c_win', `l_win' arguments should
         be positive numbers.")
  }
  if ( (abs(w_points - round(w_points)) > .Machine$double.eps ^ 0.5)
       || (abs(c_win - round(c_win)) > .Machine$double.eps ^ 0.5)
       || (abs(l_win - round(l_win)) > .Machine$double.eps ^ 0.5)){
    warning("The input values  for `w_points', `c_win', and  `l_win' should be positive integers.
            If either of them is a positive real number then the integer part of the given number
            is used to obtain the result.")
  }
  lg <- length(xd)
  w_points <- as.integer(w_points)
  c_win <- min(lg, c_win)
  c_win <- as.integer(c_win)
  l_win <- as.integer(l_win)
  if (lg <= c_win) {
    u <- L_2_norm(x = xd, thr_const = thr_con, points = w_points, rescale = rsc, gs = NULL)
    return(u)
  } else if ((lg > c_win) & (lg <= l_win)) {
    u <- L_2_norm(x = xd, thr_const = thr_con, points = w_points, rescale = rsc, gs = c_win)
    return(u)
  }
  else {
    bs <- Binary_seq(xd, grid.size = c_win)[[1]]
    cs <- cumul_sum_matrix(bs)
    t <- thr_con * sqrt(log(lg))
    K <- ceiling(lg / c_win)
    tsm <- list()
    u <- list()
    ufin <- numeric()
    uaddition <- numeric(0)
    tsm[[1]] <- xd[1:c_win]
    ufin <- L_2_norm(tsm[[1]], thr_fin = t, points = w_points, rescale = rsc)[[1]]
    if (CUSUM_at_point(x = bs,y = cs, s = max(1, c_win - (10 * w_points) + 1), e = min( (c_win + (10 * w_points)), lg), b = c_win) > t){
      uaddition <- c_win 
    }
    ufin <- c(ufin, uaddition)
    i <- 2
    while (i < K) {
      tsm[[i]] <- xd[( (i - 1) * c_win + 1):(i * c_win)]
      u[[i]] <- L_2_norm(x = tsm[[i]], thr_fin = t, points = w_points, rescale = rsc)[[1]] + (i - 1) * c_win
      if (CUSUM_at_point(x = bs,y = cs, s = max(1, i * c_win - (10 * w_points) + 1), e = min(i * c_win + (10 * w_points), lg), b = i * c_win) > t){
        uaddition <- i*c_win
      }
      ufin <- c(ufin, u[[i]], uaddition)
      i <- i + 1
    }
    tsm[[K]] <- xd[( (K - 1) * c_win + 1):lg]
    u[[K]] <- L_2_norm(tsm[[K]], thr_fin = t, points = w_points, rescale = rsc)[[1]] + (K - 1) * c_win
    ufinl <- c(ufin, u[[K]])
    return(list(cpt = sort(unique(ufinl)), bernoulli_sequence  = bs, cumul_sum = cs))
  }
}


#=================================================================
# Main functions for the Information Criterion based decision rule
#=================================================================

## The solution path obtained when the L_inf aproach is used in the overestimation step
sol_path_L_inf <- function(x, thr_ic = NULL, points = 10, resc = FALSE) {
  if (length(thr_ic) == 0){thr_ic = 0.7
  if (resc == TRUE){thr_ic = 1.6}}
  if (!(is.numeric(x))){
    stop("The input in `x' should be a numeric vector containing the data for
         which the solution path will be calculated.")
  }
  if ( (thr_ic <= 0) || (points <= 0)){
    stop("The threshold constant as well as the `points' argument that represents the
         magnitude of the expansion for the intervals should be positive numbers.")
  }
  if (abs(points - round(points)) > .Machine$double.eps ^ 0.5){
    warning("The input for `points' should be a positive integer. If it is a positive real
            number then the integer part of the given number is used as the value of `points'.")
  }
  lx_ic <- length(x)
  points <- as.integer(points)
  result_lower <- window_L_inf_norm(x, thr_con = thr_ic, w_points = points, rsc = resc)
  cpt_lower <- result_lower[[1]]
  lcpt_ic <- length(cpt_lower)
  if ((lx_ic > 1000) & (lx_ic <= 5000)){
    bs <- Binary_seq(x)[[1]]
    cs <- cumul_sum_matrix(bs)
  }
  else{bs <- result_lower[[2]]
  cs <- result_lower[[3]]
  }
  if ( (lcpt_ic == 1) | (lcpt_ic == 0)) {
    return(list(solution_path = cpt_lower, bern_sequence = bs, cum_sum = cs))
  } else {
    seb_set <- c(unique(c(1, cpt_lower)), lx_ic)
    lseb_set <- length(seb_set)
    min_C <- numeric()
    while (lseb_set >= 3) {
      Rs <- CUSUM_at_point(result_lower[[2]], y = result_lower[[3]], seb_set[1:(lseb_set - 2)], seb_set[3:(lseb_set)],
                             seb_set[2:(lseb_set - 1)], rescale = resc)
      min_Rs <- seb_set[2:(lseb_set - 1)][which.min(Rs)]
      min_C <- c(min_C, min_Rs)
      seb_set <- seb_set[-which(seb_set == min_Rs)]
      lseb_set <- lseb_set - 1
    }
    return(list(solution_path = min_C[length(min_C):1], bern_sequence = bs, cum_sum = cs))
  }
}

## The solution path obtained when the L_2 aproach is used in the overestimation step
sol_path_L_2 <- function(x, thr_ic = NULL, points = 10, resc = FALSE) {
  if (length(thr_ic) == 0){thr_ic = 0.45
  if (resc == TRUE){thr_ic = 0.8}}
  if (!(is.numeric(x))){
    stop("The input in `x' should be a numeric vector containing the data for
         which the solution path will be calculated.")
  }
  if ( (thr_ic <= 0) || (points <= 0)){
    stop("The threshold constant as well as the `points' argument that represents the
         magnitude of the expansion for the intervals should be positive numbers.")
  }
  if (abs(points - round(points)) > .Machine$double.eps ^ 0.5){
    warning("The input for `points' should be a positive integer. If it is a positive real
            number then the integer part of the given number is used as the value of `points'.")
  }
  lx_ic <- length(x)
  points <- as.integer(points)
  result_lower <- window_L_2_norm(x, thr_con = thr_ic, w_points = points, rsc = resc)
  cpt_lower <- result_lower[[1]]
  lcpt_ic <- length(cpt_lower)
  if ((lx_ic > 1000) & (lx_ic < 5000)){
    bs <- Binary_seq(x)[[1]]
    cs <- cumul_sum_matrix(bs)
  }
  else{bs <- result_lower[[2]]
  cs <- result_lower[[3]]
  }
  if ( (lcpt_ic == 1) | (lcpt_ic == 0)) {
    return(list(solution_path = cpt_lower, bern_sequence = bs, cum_sum = cs))
  } else {
    seb_set <- c(unique(c(1, cpt_lower)), lx_ic)
    lseb_set <- length(seb_set)
    min_C <- numeric()
    while (lseb_set >= 3) {
      Rs <- CUSUM_at_point(result_lower[[2]], y = result_lower[[3]], seb_set[1:(lseb_set - 2)], seb_set[3:(lseb_set)],
                             seb_set[2:(lseb_set - 1)], rescale = resc)
      min_Rs <- seb_set[2:(lseb_set - 1)][which.min(Rs)]
      min_C <- c(min_C, min_Rs)
      seb_set <- seb_set[-which(seb_set == min_Rs)]
      lseb_set <- lseb_set - 1
    }
    return(list(solution_path = min_C[length(min_C):1], bern_sequence = bs, cum_sum = cs))
  }
}


## The result of the non-parametric information criterion approach when the L_inf
## norm is employed for the first, overestimation based step.
model_selection_L_inf <- function(x, th_const = NULL, Kmax = 200, points = 10, rescale = FALSE, zeta = (log(length(x)))^(2.1)/2) {
  if (length(th_const) == 0){th_const = 0.7
  if (rescale == TRUE){th_const = 1.6}}
  print(th_const)
  if (!(is.numeric(x))){
    stop("The input in `x' should be a numeric vector containing the data in
         which you want to look for change-points.")
  }
  if ( (th_const <= 0) || (points <= 0)){
    stop("The threshold constant as well as the `points' argument that represents the
         magnitude of the expansion for the intervals should be positive numbers.")
  }
  if (abs(points - round(points)) > .Machine$double.eps ^ 0.5){
    warning("The input for `points' should be a positive integer. If it is a positive real
            number then the integer part of the given number is used as the value of `points'.")
  }
  result <- list()
  lx <- length(x)
  if (Kmax == 0 || lx <= 3) {
    result$cpt_ic <- NA
    result$no_cpt_ic <- 0
    result$ic_curve <- NA
    if (Kmax == 0) {
      stop("No change-points found, choose larger Kmax")
    } else {
      stop("Sample size is too small")
    }
  } else {
    s_p <- sol_path_L_inf(x, thr_ic = th_const, points = points, resc = rescale)
    cpt_cand <- s_p[[1]]
    b_seq <- s_p[[2]]
    cms_matrix <- s_p[[3]]
    if (length(cpt_cand) == 0){
      result$cpt_ic <- 0
      result$no_cpt_ic <- 0
      result$ic_curve <- NA
    }
    if (length(cpt_cand) > min(Kmax, lx - 2)){
      cpt_cand <- cpt_cand[1:min(Kmax, lx - 2)]
    }
    len_cpt <- length(cpt_cand)
    result$sol_path <- cpt_cand
    result$ic_curve <- rep(0, len_cpt + 1)
    if (len_cpt){
      for (i in len_cpt:1) {
        min_log_lik <- - loglik_values(x, b1 = b_seq, y = cms_matrix, cpt = cpt_cand[1:i])
        result$ic_curve[i+1] <- min_log_lik + i*zeta
      }
    }
    result$ic_curve[1] <- - loglik_values(x, b1 = b_seq, y = cms_matrix, cpt = NULL)
    result$cpt_ic <- list()
    tmp <- stats::quantile(which.min(result$ic_curve), 0.5, type = 3)
    if (tmp == 1) {
      result$cpt_ic <- NA
      result$no_cpt_ic <- as.integer(0)
    } else {
      result$cpt_ic <- sort(cpt_cand[1:(tmp - 1)])
      result$no_cpt_ic <- as.integer(tmp - 1)
    }
  }
  return(result)
}

## The result of the non-parametric information criterion approach when the L_2
## norm is employed for the first, overestimation based step.
model_selection_L_2 <- function(x, th_const = NULL, Kmax = 200, points = 10, rescale = FALSE, zeta = (log(length(x)))^(2.1)/2) {
  if (length(th_const) == 0){th_const = 0.45
  if (rescale == TRUE){th_const = 0.8}}
  if (!(is.numeric(x))){
    stop("The input in `x' should be a numeric vector containing the data in
         which you want to look for change-points.")
  }
  if ( (th_const <= 0) || (points <= 0)){
    stop("The threshold constant as well as the `points' argument that represents the
         magnitude of the expansion for the intervals should be positive numbers.")
  }
  if (abs(points - round(points)) > .Machine$double.eps ^ 0.5){
    warning("The input for `points' should be a positive integer. If it is a positive real
            number then the integer part of the given number is used as the value of `points'.")
  }
  result <- list()
  lx <- length(x)
  if (Kmax == 0 || lx <= 3) {
    result$cpt_ic <- NA
    result$no_cpt_ic <- 0
    result$ic_curve <- NA
    if (Kmax == 0) {
      stop("No change-points found, choose larger Kmax")
    } else {
      stop("Sample size is too small")
    }
  } else {
    s_p <- sol_path_L_2(x, thr_ic = th_const, points = points, resc = rescale)
    cpt_cand <- s_p[[1]]
    b_seq <- s_p[[2]]
    cms_matrix <- s_p[[3]]
    if (length(cpt_cand) == 0){
      result$cpt_ic <- 0
      result$no_cpt_ic <- 0
      result$ic_curve <- NA
    }
    if (length(cpt_cand) > min(Kmax, lx - 2)){
      cpt_cand <- cpt_cand[1:min(Kmax, lx - 2)]
    }
    len_cpt <- length(cpt_cand)
    result$sol_path <- cpt_cand
    result$ic_curve <- rep(0, len_cpt + 1)
    if (len_cpt){
      for (i in len_cpt:1) {
        min_log_lik <- - loglik_values(x, b1 = b_seq, y = cms_matrix, cpt = cpt_cand[1:i])
        result$ic_curve[i+1] <- min_log_lik + i*zeta
      }
    }
    result$ic_curve[1] <- - loglik_values(x, b1 = b_seq, y = cms_matrix, cpt = NULL)
    result$cpt_ic <- list()
    tmp <- stats::quantile(which.min(result$ic_curve), 0.5, type = 3)
    if (tmp == 1) {
      result$cpt_ic <- NA
      result$no_cpt_ic <- as.integer(0)
    } else {
      result$cpt_ic <- sort(cpt_cand[1:(tmp - 1)])
      result$no_cpt_ic <- as.integer(tmp - 1)
    }
  }
  return(result)
}

#============================================================================