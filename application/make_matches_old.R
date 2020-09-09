require(Rcpp)
require(RcppArmadillo)
require(igraph)
require(magrittr)
source('network_flame_sims_old.R')
source('FLAME_bit 2.R')
source('ATE.R')
sourceCpp('subgraph_enumerate_old.cpp')
require(purrr)
library(magrittr)
library(dplyr)
# library(splitstackshape)

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
  train <- train[, -ncol(train)]
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
  
  treated_holdout_inds <- sample(1:n_treated, round(percent_holdout * n_treated))
  
  tmp_treated_holdout <- tmp_treated[treated_holdout_inds, ]
  tmp_treated_train <- tmp_treated[-treated_holdout_inds, ]
  
  tmp_control <- filter(tmp, treated == 0)
  n_control <- nrow(tmp_control)
  
  control_holdout_inds <- sample(1:n_control, round(percent_holdout * n_control))
  
  tmp_control_holdout <- tmp_control[control_holdout_inds, ]
  tmp_control_train <- tmp_control[-control_holdout_inds, ]
  
  tmp_train <- rbind(tmp_treated_train, tmp_control_train)
  tmp_holdout <- rbind(tmp_treated_holdout, tmp_control_holdout)
  
  # If don't have both treatments and both outcomes for 
  # both train and holdout, split again
  if (length(unique(tmp_train$treated)) == 1 | 
      length(unique(tmp_train$outcome)) == 1 | 
      length(unique(tmp_holdout$treated)) == 1 |
      length(unique(tmp_holdout$outcome)) == 1) {
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

# all_dat <- load('application.RData')

## Load in data, name the adjacency matrix A, outcome Y, treatment Z, categorical covariates X
# Navigate to where the RProject is
#setwd('/Users/vittorioorlandi/Desktop/Network FLAME/Network-FLAME/')
setwd("/Users/musaidawan/Dropbox/Duke/Projects/Data for Network paper/Network-FLAME/")
village_codes <- setdiff(c(1:77), c(13,22))
net_flame_output <- matrix(NA, nrow = 77, ncol = 3)
# adjency_out <-  read.xlsx("/Users/musaidawan/Dropbox/Duke/Projects/Data for Network paper/Network-FLAME/flame_output/adjency_stats_sum_g_eq_3.xlsx",
#                          header = FALSE) %>%
#   as.matrix()
#high_d_villages <- adjency_out[which(adjency_out[,2] > 5 | adjency_out[,3] > 18)]
# adjency_out[ adjency_out[, 2] > 6 | adjency_out[, 3] > 18, ]
village_codes <- setdiff(c(9:77), c(13,22))
for (qs in c(1,2,3)) {  # c(1,2,3)
  for (val in   village_codes) { # village_codes
    if(val == 9 & qs == 3) {next}
    if(qs == 1 &  val %in% c(1,3,4,11,12,14,17,18,21,27,28,29,31,35,43,44,45,46,48, 49, 50, 55)  ) {next}
    if(qs == 2 &  val %in% c(1,4,5,6,7,8,11,40,41,42,43,44,49,50,52,54,55,56,53,65)  ) {next}
    if(qs == 3 &  val %in% c(8,9,11,40,52)  ) {next}
    #if ((val %in% high_d_villages)) {
    val <- 19
    qs <- 1
    print(paste("vilage = ", val))
    print(paste("question  = ", qs))
    # if (val %% 5 == 0) {
    #   print(val)
    # }
    #A <-
    #  read.csv(paste('./Data/Adjency/adj_andRelationships_vilno_',val,'.csv', sep = ""),
    #           header = FALSE) %>%
    #  as.matrix()

    A1 <-
      read.csv(paste('./Data/Adjency/adj_visitgo_vilno_',val,'.csv', sep = ""),
               header = FALSE) %>%
      as.matrix()

    A2 <-
      read.csv(paste('./Data/Adjency/adj_visitcome_vilno_',val,'.csv', sep = ""),
               header = FALSE) %>%
      as.matrix()

    A3 <-
      read.csv(paste('./Data/Adjency/adj_rel_vilno_',val,'.csv', sep = ""),
               header = FALSE) %>%
      as.matrix()

    A4 <-
      read.csv(paste('./Data/Adjency/adj_nonrel_vilno_',val,'.csv', sep = ""),
               header = FALSE) %>%
      as.matrix()

    #adj_visitgo_vilno_
    #adj_visitcome_vilno_
    #adj_rel_vilno_
    #adj_nonrel_vilno_

    # union of three adjency matrix
    
    A <- ((A1 + A2 + A3 + A4) >= 3) * 1
    rm(A1, A2, A3, A4)
    
    row_sum <- rowSums(A)
    print(paste('average degree is:',mean(row_sum)))
    print(paste('max degree is:',max(row_sum)))
    
    # index of units with degree > threshold
    index_max_connections <- which(row_sum > 15)
    
    demographics <- read.csv(paste('./Data/characteristics_',
                                   qs,'/village_',val,'.csv',sep =""))
    
    # row index of untreated units
    #untreated <- which(Z == 0)
    untreated <- demographics$adjmatrix_key[demographics$Z == 0]
    
    # drop units with degree connection > 15      
    demographics <- demographics[!demographics$adjmatrix_key %in% index_max_connections,]
    units_with_treatment_info <- demographics$adjmatrix_key
    #A <- A[units_with_treatment_info, units_with_treatment_info]
    Y <- demographics$Y
    Z <- demographics$Z
    X <- demographics[, which(!colnames(demographics) %in% c('adjmatrix_key', 'Y', 'Z'))]
    X <- X[,c("age","telugu")]
    rm(demographics)

    # To drop control -- control edges
    # A[untreated, untreated] <- 0
    # To drop any edges involving a control individual
    A[untreated, untreated] <- 0
    #A[, untreated] <- 0
    n <- dim(A)[1]

    # drop units with degree >threshold
    #index_max_connections <- which(row_sum > 15)
    #if (length(index_max_connections) > 0) {
    #  A <- A[-index_max_connections,]
    #  A <- A[,-index_max_connections]
    #}
        # Brute force symmetry test because isSymmetric.matrix(A) outputs FALSE
    # for (i in 1:n) {
    #   for (j in 1:n) {
    #     if (A[i, j] != A[j, i]) {
    #       print(c(i, j))
    #     }
    #   }
    # }

    # Gives graph
    G <- graph_from_adjacency_matrix(A, mode = 'undirected')

    print('Started enumerating')
    if (val %in% high_d_villages) {
      threshold <- 4
      print(paste("max threshold = ", threshold) )
      # Enumerates all possible subgraphs and puts into dataframe
      all_subgraphs <- get_neighb_subgraphs(A, units_with_treatment_info, threshold)
      
    } else {
      threshold <- 4
      print(paste("max threshold = ", threshold) )
      # Enumerates all possible subgraphs and puts into dataframe
      all_subgraphs <- get_neighb_subgraphs(A, units_with_treatment_info, threshold)
    }
    print('Finished enumerating; started classifying')
    all_features = gen_all_features(G, all_subgraphs)
    rm(all_subgraphs)
    dta = gen_data(all_features[[1]])
    saveRDS(all_features, file = paste("./flame_output/all_features/all_features_",val,'_question_',qs,".rds",sep = ""))
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
    }

    # sapply(tmp, class)

    # cols with missing values: // should be no missing values, data pre-processing
    # lapply(tmp, function(x) sum(is.na(x)) == 1)

    # If column missing < 5% of values, impute by median
    # Else, drop column
    which_missing <- which(vapply(tmp, function(x) sum(is.na(x)) > 0, 
                                  logical(1)))
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

    
    # divide tmp into training-holdout data
    split_data <- split_test_train(tmp)
    tmp_train <- split_data$train
    tmp_holdout <- split_data$holdout

    # save dta
    saveRDS(tmp_train, file = paste("./flame_output/dta/tmp_train_",val,'_question_',qs,".rds",sep = ""))
    
    
    # FLAME
    print('FLAME on!')
    flame_out <- FLAME_bit(tmp_train, tmp_holdout, 
                           A = A, network_lik_weight = 0, iterate_FLAME = TRUE)
    print('FLAME off...')
    ATE_out <- ATE(flame_out)
    print(paste("ATE = ", ATE_out))
    net_flame_output[val, qs] <- ATE_out

    df_flame_out <- flame_out$matched_data

    # write FLAME output to disk
    # Save an object to a file
    saveRDS(flame_out, file = paste("./flame_output/village_",val,'_question_',qs,".rds",sep = ""))
    # Restore the object
    #readRDS(file = "my_data.rds")
    # save(tmp_train, tmp_holdout, file = 'tmp.RData')
    #}
  }
}
