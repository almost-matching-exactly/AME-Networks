######### IMPORTANT: CHANGE TO YOUR DIRECTORY BEFORE RUNNING ############
setwd('/Users/mm567/Dropbox/Duke/projects/learning_interference/networks/Up to Date Code/application/')
set.seed(42069)

require(Rcpp)
require(RcppArmadillo)
require(igraph)
require(magrittr)
require(purrr)
library(magrittr)
library(dplyr)

source('network_flame_sims_old.R')
source('FLAME_bit 2.R')
source('ATE.R')
sourceCpp('subgraph_enumerate_old.cpp')


bin_counts = function(dta, probs = seq(0, 1, length.out = 10)){
  data.frame(apply(dta, 2, function(x){
    qs = unique(quantile(x, probs))
    if(length(unique(qs)) > 1){
      xc = as.character(cut(x, breaks=qs, include.lowest = T))
      xc[x==0] = '0'
      return(xc)
    }else{
      return(as.character(x))
    }
  }),
  stringsAsFactors = T)
}

is_degenerate <- function(train, holdout) {
  train <- train[, !(colnames(train) %in% 'matched')]
  all_dat <- 
    rbind(train, holdout) %>%
    mutate(data_type = c(rep(0, nrow(train)),
                         rep(1, nrow(holdout))))
  treatment_by_outcome <- 
    all_dat %>% 
    group_by(outcome, data_type) %>%
    summarise(n_treated = length(unique(treated))) %>%
    pull(n_treated)
  
  outcome_by_treated <- 
    all_dat %>% 
    group_by(treated, data_type) %>%
    summarise(n_outcomes = length(unique(outcome))) %>%
    pull(n_outcomes)
  
  if (any(c(treatment_by_outcome, outcome_by_treated) == 1)) {
    return(TRUE)
  }
  return(FALSE)
}

split_test_train <- function(tmp, percent_holdout = 0.1) {
  tmp_treated <- filter(tmp, treated == 1)
  n_treated <- nrow(tmp_treated)
  
  treated_holdout_inds <- sample(1:n_treated, 
                                 round(percent_holdout * n_treated))
  
  tmp_treated_holdout <- tmp_treated[treated_holdout_inds, ]
  tmp_treated_train <- tmp_treated[-treated_holdout_inds, ]
  
  tmp_control <- filter(tmp, treated == 0)
  n_control <- nrow(tmp_control)
  
  control_holdout_inds <- sample(1:n_control, 
                                 round(percent_holdout * n_control))
  
  tmp_control_holdout <- tmp_control[control_holdout_inds, ]
  tmp_control_train <- tmp_control[-control_holdout_inds, ]
  
  tmp_train <- rbind(tmp_treated_train, tmp_control_train)
  tmp_holdout <- rbind(tmp_treated_holdout, tmp_control_holdout)
  
  # If don't have both treatments and both outcomes for 
  # both train and holdout, split again
  if (is_degenerate(tmp_train, tmp_holdout)) {
    return(split_test_train(tmp, percent_holdout))
  }
  return(list(train = tmp_train,
              holdout = tmp_holdout))
}

my_combn <- function(x, m) {
  if (length(x) == 1) {
    return(list(x))
  }
  return(combn(as.integer(x), m, simplify = FALSE))
}

village_codes <- setdiff(c(1:77), c(13,22))
net_flame_output <- matrix(NA, nrow = 77, ncol = 3)

# Change t 1, 2, 3 and rerun to generate respective outputs
qs <- 1
for (val in village_codes) { 
  print(sprintf('Village %d Question %d', val, qs))
  if (val == 40) {
    next
  }
  
  A1 <-
    read.csv(paste('./Data/Adjency/adj_visitgo_vilno_',
                   val,'.csv', sep = ""),
             header = FALSE) %>%
    as.matrix()
  
  A2 <-
    read.csv(paste('./Data/Adjency/adj_visitcome_vilno_',
                   val,'.csv', sep = ""),
             header = FALSE) %>%
    as.matrix()
  
  A3 <-
    read.csv(paste('./Data/Adjency/adj_rel_vilno_',
                   val,'.csv', sep = ""),
             header = FALSE) %>%
    as.matrix()
  
  A4 <-
    read.csv(paste('./Data/Adjency/adj_nonrel_vilno_',
                   val,'.csv', sep = ""),
             header = FALSE) %>%
    as.matrix()
  
  # union of three adjency matrix
  A <- ((A1 + A2 + A3 + A4) >= 3) * 1
  rm(A1, A2, A3, A4)
  
  row_sum <- rowSums(A)
  
  # index of units with degree > threshold
  index_max_connections <- which(row_sum > 15)
  
  demographics <- read.csv(paste('./Data/characteristics_',qs,
                                 '/village_',val,'.csv',sep =""))
  
  # row index of untreated units
  #untreated <- which(Z == 0)
  untreated <- demographics$adjmatrix_key[demographics$Z == 0]
  
  # drop units with degree connection > 15      
  demographics <- demographics[!demographics$adjmatrix_key 
                               %in% index_max_connections,]
  
  units_with_treatment_info <- demographics$adjmatrix_key
  #A <- A[units_with_treatment_info, units_with_treatment_info]
  Y <- demographics$Y
  Z <- demographics$Z
  X <- demographics[, which(!colnames(demographics) 
                            %in% c('adjmatrix_key', 'Y', 'Z'))]
  X <- X[,c("age","telugu")]
  rm(demographics)
  
  n <- dim(A)[1]
  
  # Gives graph
  G <- graph_from_adjacency_matrix(A, mode = 'undirected')
  G <- set_vertex_attr(G, 'color', value = Z)#, index = units_with_treatment_info)
  
  # print('Started enumerating')
  if (max(row_sum) > 15 | mean(row_sum) > 8) {
    threshold <- 2
    # print(paste("max threshold = ", threshold) )
    # Enumerates all possible subgraphs and puts into dataframe
    all_subgraphs <- get_neighb_subgraphs(A, units_with_treatment_info,
                                          threshold)
    
  } else {
    threshold <- 4
    # print(paste("max threshold = ", threshold) )
    # Enumerates all possible subgraphs and puts into dataframe
    all_subgraphs <- get_neighb_subgraphs(A, units_with_treatment_info, 
                                          threshold)
  }
  # print('Finished enumerating; started classifying')
  all_features = gen_all_features(G, all_subgraphs)
  rm(all_subgraphs)
  dta = gen_data(all_features[[1]])
  
  saveRDS(all_features, 
          file = paste("./application_output/features/all_features_",
                       val,'_question_',qs,".rds",sep = ""))
  rm(all_features)
  # Add covariate information
  ## Check that order of covariates X is same as order of subgraph counts
  dta = bin_counts(dta, seq(0, 1, length.out = 10))  #bin_counts(dta)
  
  dta <- cbind(dta, X)
  
  # Convert everything to factor
  dta = data.frame(sapply(dta, factor), stringsAsFactors = T)
  
  # Adds outcome and treatment
  dta$outcome = as.numeric(Y) # Y should be numeric
  dta$treated = factor(Z) # Z should be binary vector
  
  # drop cols with no variation
  # should be done in FLAME, but just in case.
  drop_these <- which(lapply(dta, function(x) length(unique(x))) == 1)
  if (length(drop_these) > 0) {
    tmp <- dta[, -drop_these]
  } else {
    tmp <- dta
  }
  
  # sapply(tmp, class)
  
  # If column missing < 5% of values, impute by median
  # Else, drop column
  which_missing <- which(vapply(tmp, function(x) 
    sum(is.na(x)) > 0, logical(1)))
  if (length(which_missing) > 0) {
    to_drop <- NULL
    for (i in 1:length(which_missing)) {
      missing_col <- which_missing[i]
      missing_vals <- which(is.na(tmp[, missing_col]))
      if (length(missing_vals) / nrow(tmp) < 0.05) { # Under 5% missing
        # Impute by median
        tmp[missing_vals, missing_col] <-
          as.factor(median(as.numeric(tmp[, missing_col]), na.rm = TRUE))
      }
      else {
        to_drop <- c(to_drop, missing_col)
      }
    }
    tmp <- tmp[, -to_drop]
  }
  
  # FLAME
  print('FLAME on!')
  
  flame_out <- FLAME_bit(tmp, tmp, A = A, 
                         network_lik_weight = 0, iterate_FLAME = TRUE)
  
  # print('FLAME off...')
  ATE_out <- ATE(flame_out)
  print(paste("ATE = ", ATE_out))
  
  net_flame_output[val, qs] <- ATE_out
  
  df_flame_out <- flame_out$matched_data
  
  # write FLAME output to disk
  # Save an object to a file
  saveRDS(flame_out, file = paste("./application_output/matches/village_",
                                  val,'_question_',qs,".rds",sep = ""))
}
