# Outcome Weighting -------------------------------------------------------
LUE <- function(weights, Y) {
  return(sum(weights * Y))
}


# All Eigenvectors --------------------------------------------------------
weighted_l2 = function(x, y, W){
  sqrt(t(x - y) %*% W %*% (x - y))
}

all_eigs <- function(A, Z, Y) {
  n <- length(Z)
  eig <- eigen(A, symmetric = TRUE)
  total_treatment_effect <- 0
  ########## Might make more sense to weigh the evecs by their inverse evals, no...? 
  weights <- diag(1 / (1:ncol(eig$vectors))) # Check if eig$vectors = n 
  for (i in 1:n) {
    opposite_treatment <- (Z != Z[i])
    mi <- sapply(which(opposite_treatment),
                 function(x) {
                   weighted_l2(eig$vectors[i, ], eig$vectors[x, ], weights)})
    
    mi <- which.min(mi)
    
    Ymie = Y[opposite_treatment][mi]
    total_treatment_effect <- total_treatment_effect +
      (Y[i] - Ymie) * Z[i] + (Ymie - Y[i]) * (1 - Z[i])
  }
  return(total_treatment_effect / n)
}

# First Eigenvector -------------------------------------------------------
first_eig <- function(A, Z, Y) {
  n <- length(Z)
  total_treatment_effect <- 0
  eig <- eigen(A, symmetric = TRUE)
  ev <- eig$vectors[, order(eig$values)[1]] ## Should be unnecessary to sort
  for (i in 1:n) {
    opposite_treatment <- (Z != Z[i])
    Ymie = Y[opposite_treatment][which.min(abs(ev[i] - ev[opposite_treatment]))]
    total_treatment_effect <- total_treatment_effect +
      (Y[i] - Ymie) * Z[i] + (Ymie - Y[i]) * (1 - Z[i])
  }
  return(total_treatment_effect / n)
}

# True Estimator ----------------------------------------------------------
true_est <- function(Y, Z, f) {
  n <- length(Z)
  total_treatment_effect <- 0
  for (i in 1:n) {
    opposite_treatment <- (Z != Z[i])
    Ymif <- Y[opposite_treatment][which.min(abs(f[i] - f[opposite_treatment]))]
    total_treatment_effect <- total_treatment_effect +
      (Y[i] - Ymif) * Z[i] + (Ymif - Y[i]) * (1 - Z[i])
  }
  return(total_treatment_effect / n)
}

# Naive Stratified Degree Estimator ---------------------------------------

n_fun <- function(z, d, Z, treated_degree) {
  sum(Z == z & treated_degree == d)
}

C_strat_naive <- function(d, Z, treated_degree) {
  n_0d <- n_fun(0, d, Z, treated_degree)
  n_1d <- n_fun(1, d, Z, treated_degree)
  if (n_0d * n_1d == 0) {
    return(0)
  }
  return(1 / (1 / n_0d + 1 / n_1d))
}
  
strat_naive <- function(A, Z, Y) {
  ## Implements stratified naive estimator as given in 
  ## Proposition 6.5 of Sussman & Airoldi 2017
  d <- colSums(A)
  treated_degree <- colSums(A * Z)
  n <- dim(A)[1]
  all_C <- vapply(1:n, function(i) {C_strat_naive(d[i], Z, treated_degree)}, 
                  FUN.VALUE = numeric(1))
  C_weight <- all_C / sum(all_C)
  weights <- sapply(1:n, function(i) {
    a <- C_weight[i]
    b <- (2 * Z[i] - 1) / n_fun(Z[i], treated_degree[i], Z, treated_degree)
    a * b
  })
  return(LUE(weights, Y))
}

# (Independent) SANIA -----------------------------------------------------

C_sum_SANIA <- function(n_neighbs, Z, p) {
  # Need to code this more efficiently 
  
  n <- length(Z)
  C_sum <- vector(mode = 'numeric', length = n)
  sq_p <- Z * p ^ 2 + (1 - Z) * (1 - p) ^ 2
  for (i in 1:n) {
    tmp <- 0
    for (d in 0:n_neighbs[i]) {
      tmp <- tmp + choose(n_neighbs[i], d) ^ 2 * p ^ d * (1 - p) ^ (n_neighbs[i] - d)
    }
    C_sum[i] <- tmp * sq_p[i]
  }
  return(C_sum)
}

ind_SANIA <- function(G, Z, Y, p) {
  ## Implements the MIV LUE in the SANIA setting where priors are independent. 
  ## Specified in Theorem 6.2 of Sussman and Airoldi 2017.
  ## Currently assumes treatment assigned Bernoulli(p)
  
  n <- length(Z)
  n_neighbs <- ego_size(G) - 1
  C_sum <- C_sum_SANIA(n_neighbs, Z, p)
  weights <- (Z * p - (1 - Z) * (1 - p)) / (n * C_sum)
  return(LUE(weights, Y))
}

get_Cid <- function(P_zd, p, max_deg, n_neighbs, Z) {
  n <- nrow(A)
  Cid <- matrix(0, nrow = max_deg + 1, ncol = n)
  for (i in 1:n) {
    for (d in 0:n_neighbs[i]) {
      if (P_zd[d + 1, Z[i] + 1] == 0) {
        Cid[d + 1, n] <- 0
      }
      else {
      Cid[d + 1, n] <- 
        choose(n_neighbs[i], d) * p ^ d  * (1 - p) ^ (n_neighbs[i] - d) / 
        (P_zd[d + 1, Z[i] + 1]) ^ 2
      }
    }
  }
  return(Cid)
}

naive_SANIA <- function(G, A, Z, Y, p) {
  # AZ <- A
  # AZ[which(Z == 0), ] <- 0
  # AZ[, which(Z == 0)] <- 0
  treated_degree <- rowSums(A) # Assuming treated A
  # treated_degree <- rowSums(AZ) # Assuming treated A
  n <- nrow(A)
  n_neighbs <- ego_size(G) - 1
  max_deg <- max(treated_degree)
  P_zd <- matrix(, nrow = max_deg + 1, ncol = 2)
  for (d in 0:max_deg) {
    for (z in 0:1) {
      same_deg <- which(treated_degree == d)
      P_zd[d + 1, z + 1] = sum(Z[same_deg] == z) / n
    }
  }
  Cid <- get_Cid(P_zd, p, max_deg, n_neighbs, Z)
  Cid_sum <- colSums(Cid) 
  weights <- vector(mode = 'numeric', length = n)
  for (i in 1:length(weights)) {
    if (P_zd[treated_degree[i] + 1, Z[i] + 1] == 0 | 
        Cid_sum[i] == 0) {
      weights[i] <- 0
    }
    else {
      weights[i] <- 
        (2 * Z[i] - 1) * Cid[treated_degree[i], i] / 
        (n * Cid_sum[i] * P_zd[treated_degree[i], Z[i] + 1])
    }
  }
  return(LUE(weights, Y))
}

# SANASIA -----------------------------------------------------------------
### Incomplete 
sigma <- function(sig_sq_a, sig_sq_b, sig_sq_c, A, Z) {
  ## For computing SANASIA estimator from Airoldi & Sussman 2017
  d_z <- colSums(A * Z)
  return(sig_sq_a + sig_sq_b * Z + sig_sq_c * d_z * sum(d_z))
}

SANASIA <- function(A, Z, Y, treatment_prob, treatment_num, sig_sq_a, sig_sq_b, sig_sq_c) {
  n <- dims(A)[1]
  n_treated <- sum(Z)
  if (is.null(treatment_num)) {
    p_z <- treatment_prob ^ n_treated * 
      (1 - treatment_prob) ^ (n - n_treated)
    
  }
  return(LUE(weights, Y))
}
