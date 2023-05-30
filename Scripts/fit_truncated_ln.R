require(tidyverse)
require(cubature)

# Following method of
# Pueyo S (2014) Algorithm for the maximum likelihood estimation of the parameters
# of the truncated normal and lognormal distributions 23–6 http://arxiv.org/abs/1407.6518

# User provides:
# --- patch_dist: patch size distribution as vector
# --- xmin: lower truncation value
# --- xmax: upper truncation value
# --- eta: convergence parameter; 0 < eta < 1

fit_truncated_ln <- function(patch_dist, xmin, xmax = Inf, eta = 0.1){
  
  # Filter patch size distribution to only include values >= xmin
  patch_dist <- patch_dist[patch_dist >= xmin]
  
  # Log-transform data
  y <- log(patch_dist)
  ymin <- log(xmin)
  ymax <- log(xmax)
  
  # Calculate sampling means
  y_bar <- mean(y)
  y2_bar <- mean(y^2)
  y3_bar <- mean(y^3)
  y4_bar <- mean(y^4)
  
  # Calculate h, a, b, c
  h <- y4_bar*(-y2_bar + y_bar^2) + y3_bar*(y3_bar - 2*y_bar*y2_bar) + y2_bar^3
  a <- (y4_bar - y2_bar^2) / h
  b <- (-y3_bar + y_bar*y2_bar) / h
  c <- (y2_bar - y_bar^2) / h
  
  # Starting values for alpha and psi
  alpha_j <- 1
  psi_j <- 0.1
  
  # Starting values for deltas
  delta_j_alpha <- delta_j_psi <- 1
  
  ##### Implement algorithm
  while (abs(delta_j_alpha) > eta/100000 & abs(delta_j_psi) > eta/100000) {

    # Normalization constant
    f_con <- function(u) {exp(-alpha_j*u - psi_j*u^2)}
    constant <- cubintegrate(f_con, lower = ymin, upper = ymax, method = "pcubature")$integral
    
    # Expectation of y
    f_Ey <- function(y) {y * (exp(-alpha_j*y - psi_j*y^2) / constant)}
    Ey <- cubintegrate(f_Ey, lower = ymin, upper = ymax, method = "pcubature")$integral
    
    # Expectation of y^2
    f_Ey2 <- function(y) {y^2 * (exp(-alpha_j*y - psi_j*y^2) / constant)}
    Ey2 <- cubintegrate(f_Ey2, lower = ymin, upper = ymax, method = "pcubature")$integral
    
    # Calculate deltas
    delta_j_alpha <- a*eta*(y_bar - Ey) + b*eta*(y2_bar - Ey2)
    delta_j_psi <-   b*eta*(y_bar - Ey) + c*eta*(y2_bar - Ey2)
    
    # Update of parameters
    alpha_j <- alpha_j + delta_j_alpha
    psi_j <- psi_j + delta_j_psi
  }
  
  # Return final parameters
  return( list(beta = 1 + alpha_j, psi = psi_j) )

}


