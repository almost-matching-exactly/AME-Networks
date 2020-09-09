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

drop_unmatchable <- function(data) {
  # Called before every covariate is temporarily dropped and the result analyzed
  # Drops any covariates where there are two levels of a factor, one of which is only
  #   exhibited by a single unit. Also 1-level covariates
  # Will think more about cases where there are n levels of a factor, n - 1 (k?) of which
  #   are only exhibited by a single unit
  # Makes the implicit assumption that these covariates are not useful for prediction
  # print("in here")
  dims <- dim(data)
  if (dims[1] == 2) {
    return(data)
  }
  cov_data <- data[, 1:(dims[2] - 2)]
  unmatchable <- which(vapply(cov_data, function(x) {
    (length(unique(x)) == 1) | ((length(unique(x)) == 2) & (1 %in% table(x)))},
    FUN.VALUE = logical(1)))
  n_unmatchable <- length(unmatchable)
  if (n_unmatchable > 0) {
    # print(sprintf('number unmatchable is %d', n_unmatchable))
    return(list(data = data[, -unmatchable],
                unmatchable = unmatchable))
  }
  return(list(data = data,
              unmatchable = c()))
}

aggregate_table <- function(list) {
  tab = table(as.character(list))
  tab = unclass(tab)
  name = names(tab)
  list_val = as.character(list)
  return(as.vector(tab[sapply(list_val, function(x) which(name ==  x))]))
}

# update_matched_bit takes a dataframe, a set of covariates to match on,
# the treatment indicator column and the matched indicator column.
# it returns the array indicating whether each unit is matched (the first return value),
# and a list of indices for the matched units (the second return value)

update_matched_bit <- function(data, cur_covs, covs_max_list, compute_var) {

  # browser()
  ## Tentatively going to say the following:
  # p <- ncol(data)
  # data_wo_t <- as.bigz(as.matrix(data[, 1:(p - 3)])) # -1 for outcome, -1 for matched, -1 for treatment
  data_wo_t <- as.bigz(as.matrix(data[, cur_covs+1])) # the covariates values as a matrix

  options("scipen" = 100, "digits" = 4)

  # Compute b_u
  multiplier = as.bigz(rep(0, length(cur_covs)))
  for (i in 1:length(cur_covs)) {
    multiplier[i] = pow.bigz(covs_max_list[i], i - 1)
  }
  # multiplier is in fact a bigz
  b_u = as.vector(data_wo_t %*% as.matrix(multiplier))

  # Compute b_u+
  multiplier = as.bigz(rep(0,length(cur_covs)))
  for (i in 1:length(cur_covs)) {
    multiplier[i] = pow.bigz(covs_max_list[i], i)
  }

  b_u_plus = as.vector(data_wo_t %*% as.matrix(multiplier))
  b_u_plus = add.bigz(b_u_plus, data[,'treated'])

  # Compute c_u
  c_u = aggregate_table(b_u)

  # Compute c_u+
  c_u_plus = aggregate_table(b_u_plus)

  if (compute_var) {
    match_index = mapply(function(x,y) (x != y) && (x >= 4) && (y >= 2) && (x - y >= 2), c_u, c_u_plus)
  }
  else {
    match_index = mapply(function(x,y) (x != y) && (x >= 2) && (y >= 1), c_u, c_u_plus)
  }
  index = b_u[match_index]
  return(list(match_index, index))
}

#get_CATE function takes match_index and index (b_u values)
# and return dataframe that includes
#(1) list of covariates that are used to match at level l
#(1) conditional average treatment effect (effect)
#(2) size of each matched group (size)

get_CATE_bit <- function(data, match_index, index, cur_covs, covs_max_list, column, factor_level, compute_var, num_covs) {
  if (length(index) == 0) {
    if (compute_var) {
      CATE <- setNames(data.frame(matrix(ncol = length(cur_covs)+3, nrow = 0)),
                       c(column[(cur_covs + 1)],"effect","size", "variance"))
    }
    else {
      CATE <- setNames(data.frame(matrix(ncol = length(cur_covs)+2, nrow = 0)),
                       c(column[(cur_covs + 1)],"effect","size"))
    }
  }

  else {
    # browser()
    d = data[match_index,]
    d[,'b_u'] = index
    d[,'b_u'] = unlist(lapply(d[,'b_u'], as.character))
    d[,1:num_covs] <- mapply(function(x,y) factor_level[[x]][d[,y]], 1:num_covs, 1:num_covs)

    summary = data.frame(d %>% group_by(.data$b_u,.data$treated) %>%
                           summarise(size = length(.data$outcome), mean = mean(.data$outcome), variance= var(.data$outcome)))
    summary = data.frame(summary %>% group_by(.data$b_u) %>%
                           summarize(size = sum(.data$size), treated_lst = list(.data$treated), mean_lst = list(.data$mean), var_list = list(.data$variance)))


    pos <- unlist(lapply(summary$b_u, function(x) which(d$b_u %in% x)[1]))
    CATE <- as.data.frame(d[pos, cur_covs + 1])


    CATE$effect = mapply(function(x,y) x[which(y == 1)] - x[which(y == 0)], summary$mean_lst, summary$treated_lst)
    CATE$size = summary$size

    if (compute_var) {
      CATE$variance = mapply(correct_variance, summary$var_list, summary$treated_lst)
      colnames(CATE) = c(column[(cur_covs + 1)],"effect","size", "variance")
    }
    else {
      colnames(CATE) = c(column[(cur_covs + 1)],"effect","size")
    }

    CATE <- CATE[order(CATE$effect),]
    rownames(CATE) = NULL
  }
  return(CATE)
}

correct_variance <- function(x,y) {
  if (is.null(x)) {
    return(0)
  }
  else {
    return((x[which(y == 1)]) + (x[which(y == 0)]))
  }
}

Regression_PE_bit <- function(holdout_trt, holdout_ctl) {

  # MSE for treated
  model_lm <- lm(outcome ~ ., data = holdout_trt) # fit the data to lm model
  MSE_treated <- mean((holdout_trt$outcome - model_lm$fitted.values)^2) # compute mean squared error

  # MSE for control
  model_lm <- lm(outcome ~ ., data = holdout_ctl) # fit the data to lm model
  MSE_control <- mean((holdout_ctl$outcome - model_lm$fitted.values)^2)# compute mean squared error

  return(MSE_treated + MSE_control)
}

GLMNET_PE_bit <- function(holdout_trt, holdout_ctl, lambda, alpha) {

  # MSE for treated
  y <- holdout_trt$outcome
  x <- model.matrix(~ .-1, holdout_trt[,-which(colnames(holdout_trt) == "outcome")])
  #print(apply(x, 2, function(j) sum(j != 0) ))
  #x <- x[, apply(x, 2, var, na.rm=TRUE) > 0]
  fit <- glmnet(x, y, alpha = alpha, lambda = lambda)
  predicted_value <- predict(fit, x, s = lambda)
  MSE_treated <- mean((y - predicted_value)^2) # compute mean squared error

  # MSE for control
  y <- holdout_ctl$outcome
  x <- model.matrix(~ .-1, holdout_ctl[,-which(colnames(holdout_ctl) == "outcome")])
  #x <- x[, apply(x, 2, var, na.rm=TRUE) > 0]
  fit <- glmnet(x, y, alpha = alpha, lambda = lambda)
  predicted_value <- predict(fit, x, s = lambda)
  MSE_control <- mean((y - predicted_value)^2) # compute mean squared error

  return(MSE_treated + MSE_control)
}


#match_quality function takes holdout dataset, number of total covariates,
#list of current covariates, covariate c to temporily remove from, and trafeoff
#parameter as input. The function then computes Balancing Factor and Predictive Error,
#returning Match Quality.

match_quality_bit <- function(c, data, holdout, num_covs, cur_covs, covs_max_list, tradeoff,
                              PE_function, model, ridge_reg, lasso_reg, compute_var,
                              A = NULL, network_lik_weight = 0) {


  # browser()
  # temporarily remove covariate c
  covs_to_match = cur_covs[cur_covs != c]
  covs_max_to_match = covs_max_list[-which(cur_covs == c)]

  # Calculate number of units unmatched (available)

  num_control = nrow(data[data[,'treated'] == 0,])
  num_treated = nrow(data[data[,'treated'] == 1,])

  # Number of matched units

  match_index = update_matched_bit(data, covs_to_match, covs_max_to_match, compute_var)[[1]]

  num_control_matched = nrow(data[match_index & data[,'treated'] == 0,])
  num_treated_matched = nrow(data[match_index & data[,'treated'] == 1,])

  # Compute Predictive Error

  holdout_trt <- holdout[holdout[,'treated'] == '1',-(c+1)]
  holdout_trt <- holdout_trt[,!(names(holdout_trt) %in% 'treated')]
  holdout_ctl <- holdout[holdout[,'treated'] == '0',-(c+1)]
  holdout_ctl <- holdout_ctl[,!(names(holdout_ctl) %in% 'treated')]

  ## Eliminate 0-variance predictors. Sloppy, need to fix.
  # all_same <- which(sapply(holdout_ctl, function(x) length(unique(x)) == 1))
  # if (length(all_same) > 0) {
  #   holdout_ctl <- holdout_ctl[, -all_same]
  # }
  # all_same <- which(sapply(holdout_trt, function(x) length(unique(x)) == 1))
  # if (length(all_same) > 0) {
  #   holdout_trt <- holdout_trt[, -all_same]
  # }

  if (is.null(PE_function)) {

    # default PE - ridge regression with 0.1 regularization parameter
    if (is.null(model)) {
      # browser()
      PE <- GLMNET_PE_bit(holdout_trt, holdout_ctl, lambda = 0.1, alpha = 0)
    }
    else {
      if (model == "Linear") {
        PE <- Regression_PE_bit(holdout_trt, holdout_ctl)
      }

      if (model == "Lasso") {
        if (is.null(lasso_reg)) {
          stop("Please specify lasso_reg regularization parameter.")
        }
        PE <- GLMNET_PE_bit(holdout_trt, holdout_ctl, lambda = lasso_reg, alpha = 1)
      }

      if (model == "Ridge") {
        if (is.null(ridge_reg)) {
          stop("Please specify ridge_reg regularization parameter")
        }
        PE <- GLMNET_PE_bit(holdout_trt, holdout_ctl, lambda = ridge_reg, alpha = 0)
      }
    }
  }

  else {
    # Compute PE based on user defined PE_function
    PE <- PE_function(holdout_trt$outcome, holdout_ctl$outcome, cbind(holdout_trt[,-which(colnames(holdout_trt) == "outcome")]),
                      cbind(holdout_ctl[,-which(colnames(holdout_ctl) == "outcome")]))
  }

  BF <- if_else(num_control == 0 | num_treated == 0,
                0,
                num_control_matched / num_control + num_treated_matched / num_treated)

  if (network_lik_weight == 0) {
    return(tradeoff * BF - PE)
  }

  X <- as.matrix(rbind(holdout_ctl, holdout_trt))
  X <- X[, -ncol(X)]
  p <- ncol(X)

  edges <- A[upper.tri(A)]
  edge_mat <- matrix(, nrow = length(edges), ncol = 2 * p)
  n_units <- dim(A)[1]
  all_pairs <- combn(c(1:n_units), 2, simplify = FALSE)

  for (i in seq_along(all_pairs)) {
    edge_mat[i, ] <- cbind(X[all_pairs[[i]][1], ], X[all_pairs[[i]][2], ])
  }
  edge_df <- cbind(as.data.frame(edge_mat), edges)
  # print('got here!')
  all_same <- which(sapply(edge_df, function(x) length(unique(x)) == 1))
  if (length(all_same) > 0) {
    edge_df <- edge_df[, -all_same]
  }
  logistic_model <- glm(edges ~ ., family = 'binomial', data = edge_df)
  # print('and here!')
  log_lik <- logLik(logistic_model)
  # print(log_lik)
  return(tradeoff * BF + log_lik * network_lik_weight - PE)
}


#' Bit Vectors Implementation
#'
#' \code{FLAME_bit} applies FLAME matching algorithm based on bit vectors.
#' The required arguments include (1) data and (2) holdout. The default model
#' for Match Quality is set to Ridge regression with 0.1 regularization parameter.
#'
#' @param data input data
#' @param holdout holdout training data
#' @param compute_var variance indicator (optional, default = FALSE)
#' @param tradeoff Match Quality tradeoff parameter (optional, default =
#'   0.1)
#' @param PE_function user defined function to compute predictive error
#'   (optional)
#' @param model Linear, Ridge, or Lasso (optional)
#' @param ridge_reg L2 regularization parameter if model = Ridge (optional)
#' @param lasso_reg L1 regularization parameter if model = Lasso (optional)
#' @return (1) list of covariates FLAME performs matching at each iteration, (2)
#' Sizes, conditional average treatment effects (CATEs), and variance (if compute_var = TRUE)
#' of matches at each iteration, (3) match quality at each iteration, and (4) the original
#' data with additional column *matched*, indicating the number of covariates each unit is
#' matched on. If a unit is never matched, then *matched* will be 0.
#' @examples
#' data(toy_data)
#' FLAME_bit(data = toy_data, holdout = toy_data)
#' @import dplyr
#' @import gmp
#' @import glmnet
#' @importFrom rlang .data
#' @importFrom graphics boxplot
#' @importFrom stats rbinom rnorm runif setNames
#' @importFrom stats lm var
#' @export

FLAME_bit <- function(data, holdout, tradeoff = 0.1, compute_var = FALSE, PE_function = NULL,
                      model = NULL, ridge_reg = NULL, lasso_reg = NULL,
                      A = NULL, network_lik_weight = 0, iterate_FLAME = FALSE) {

  require(stringr)
  require(gmp)
  require(glmnet)
  require(dplyr)

  num_covs <- ncol(data) - 2 # ignore treatment and outcome
  
  # Stop if covariates are not factors
  if (Reduce("|", sapply(1:num_covs, function(x) !is.factor(data[,x] )))) {
    stop("Covariates are not factor data type.")
  }
  for (val in 1:num_covs ) {
    if (  !is.factor(holdout[,val])  ) {
      stop("Covariates are not factor data type.")
    }
  }
  
#  if (Reduce("|", sapply(1:num_covs, function(x) !is.factor(data[,x] ))) |
#      Reduce("|", sapply(1:num_covs, function(x) !is.factor(holdout[,x] )))) {
#    stop("Covariates are not factor data type.")
#  }

  # Stop if treatment isn't factor
  if (!is.factor(data[,num_covs + 2]) | !is.factor(holdout[,num_covs + 2])) {
    stop("Treatment variable is not factor data type")
  }

  # Stop if outcome is not numeric
  if (!is.numeric(data[,num_covs + 1]) | !is.numeric(holdout[,num_covs + 1])) {
    stop("Outcome variable is not numeric data type")
  }
  
  unm_data <- drop_unmatchable(data)$unmatchable
  unm_hld <- drop_unmatchable(holdout)$unmatchable
  unmatchable <- union(unm_data, unm_hld)
  if (!is.null(unmatchable)){
    data <- data[, -unmatchable]
    holdout <- holdout[, -unmatchable]
  }
  num_covs <- ncol(data) - 2 # ignore treatment and outcome

  factor_level <- lapply(data[, 1:num_covs], levels)  # Levels of each factor
  covs_max_list <- sapply(factor_level, length) # Number of levels of each covariate
  
  # Sort in increasing order of number of levels
  covs_max_list <- covs_max_list[order(covs_max_list)]
  factor_level <- factor_level[names(covs_max_list)]
  
  data[, c(1:num_covs)] <- data[,names(covs_max_list)]
  colnames(data) <- c(names(covs_max_list), "outcome", "treated")
  holdout[, c(1:num_covs)] <- holdout[,names(covs_max_list)]
  colnames(holdout) <- c(names(covs_max_list), "outcome", "treated")
  
  # Add column "matched" to input data
  data$matched <- as.integer(0)
  column <- colnames(data)

  # Convert each covariate and treated into type integer
  data[,c(1:num_covs, num_covs + 2)] <- sapply(data[,c(1:num_covs, num_covs + 2)], function(x) as.integer(x))
  data$treated <- data$treated - 1 # Treatment goes to 1 and 2 for untreated and treated when made into integer
  
  # Change input data and holdout training data column name
  colnames(data) <- c(paste("x", seq(0, num_covs - 1), sep = ""), "outcome", "treated", "matched")
  colnames(holdout) <- c(paste("x", seq(0, num_covs - 1), sep = ""), "outcome", "treated")

  #Set up return objects
  covs_list <- list() # List of covariates for matching at each level
  CATE <- list() # List of dataframe that calculates conditional average treatment effect at each level
  SCORE <- list()
  matched_groups <- list()
  n_matched_groups <- 0

  # Initialize the current covariates to be all covariates and set level to 1
  cur_covs <- seq(0, num_covs - 1)
  level <- 1
  
  # Get matched units without dropping anything
  return_match <- update_matched_bit(data, cur_covs, covs_max_list, compute_var)
  match_index <- return_match[[1]]
  index <- return_match[[2]]

  # Set matched = num_covs and get those matched units
  data[match_index,'matched'] = length(cur_covs)
  return_df = data[match_index,]

  covs_list[[level]] <- column[(cur_covs + 1)]
  CATE[[level]] <- get_CATE_bit(data, match_index, index, cur_covs, covs_max_list, column, factor_level, compute_var, num_covs)

  data = data[!match_index,]

  # While there are still covariates for matching
  if (iterate_FLAME) {
    n_cov_min <- 1
  }
  else {
    n_cov_min <- Inf
  }
  my_count <- 0
  
  while (length(cur_covs) > n_cov_min && # while((length(cur_covs) > 1)...)
         (sum(data[,'treated'] == 0) > 0) &&
         (sum(data[,'treated'] == 1) > 0)) {
    
    if (is_degenerate(data, holdout)) {
      break
    }
    
    my_count <- my_count + 1
    unm_data <- drop_unmatchable(data)$unmatchable
    unm_hld <- drop_unmatchable(holdout)$unmatchable
    unmatchable <- union(unm_data, unm_hld)
    if (!is.null(unmatchable)){
     cur_covs <- cur_covs[-unmatchable]
     covs_max_list <- covs_max_list[-unmatchable]
     if (length(cur_covs) == 0) {
       break
     }
    }
    
    # print('Data variance')
    # print(apply(data, 2, function(j) c(var(j), var(j[data$treated==1]), var(j[data$treated==0]))))
    # print('Holdout variance')
    # print(apply(holdout, 2, function(j) c(var(j), var(j[holdout$treated==1]), var(j[holdout$treated==0]))))
    # browser()
    level = level + 1

    #Temporarily drop one covariate at a time to calculate Match Quality
    #Drop the covariate that returns highest Match Quality Score
    
    list_score <- unlist(lapply(cur_covs, match_quality_bit, data, holdout, num_covs, cur_covs, covs_max_list,
                                tradeoff, PE_function, model, ridge_reg, lasso_reg, compute_var,
                                A, network_lik_weight)) # Only last 2 args are new
    quality <- max(list_score)

    # randomly sample one covariate to drop
    if (length(quality) > 1) {
      drop <- sample(which(list_score == quality),1)
    }
    else {
      drop <- which(list_score == quality)
    }

    cur_covs <- cur_covs[-drop]

    if (length(cur_covs) == 0) {
      break
    }

    covs_max_list = covs_max_list[-drop]

    SCORE[[level-1]] <- quality
    covs_list[[level]] <- column[(cur_covs + 1)]

    # Update Match

    return_match = update_matched_bit(data, cur_covs, covs_max_list, compute_var)
    match_index = return_match[[1]]
    index = return_match[[2]]

    ## For ATE computation
    # if (sum(match_index) > 0) {
    #   matched_inds <- which(match_index)
    #   for (i in seq_along(unique(index))) {
    #     MG <- matched_inds[which(index == index[i])]
    #     MG_outcomes <- data$outcome[MG]
    #     MG_treatments <- data$treated[MG]
    #     MG_control <- MG_treatments == 0
    #     MG_treated <- MG_treatments == 1
    #     mean_control_outcome <- mean(MG_outcomes[MG_control]) # Avg control outcome in this MG
    #     TTT <- TTT + sum(MG_outcomes[MG_treated]) - sum(MG_treated) * mean_control_outcome
    #     n_matched_treated <- n_matched_treated + sum(MG_treated)
    #   }
    # }

    # Set matched = num_covs and get those matched units
    data[match_index,'matched'] = length(cur_covs)
    
    sum(match_index)
    # browser()
    return_df = rbind(return_df,data[match_index,])
    #browser()
    CATE[[level]] <- get_CATE_bit(data, match_index, index, cur_covs, covs_max_list, column, factor_level, compute_var, num_covs)

    # if (sum(match_index) > 0) { # Have a new matched group
    #   MG_units <- which(match_index) # Units matched this iteration
    #   MG_outcomes <- data$outcome[MG_units]
    #   MG_treatments <- data$treated[MG_units]
    #   MG_control <- MG_treatments == 0
    #   MG_treated <- MG_treatments == 1
    #   n_matched_treated <- n_matched_treated + sum(MG_treated)
    #   mean_control_outcome <- mean(MG_outcomes[MG_control]) # Avg control outcome in this MG
    #   TTT <- TTT + sum(MG_outcomes[MG_treated]) - n_matched_treated * mean_control_outcome
    # }
    # Remove matched_units
    data = data[!match_index,]
    # message(paste("number of matched units =", sum(match_index)))

  }
  if (nrow(data) != 0) {
    return_df = rbind(return_df,data)
  }
  colnames(return_df) <- column
  rownames(return_df) <- NULL
  return_df[,1:num_covs] <- mapply(function(x,y) factor_level[[x]][return_df[,y]], 1:num_covs, 1:num_covs)
  return_df$index <- 1:nrow(return_df)
  return_list = list(covs_list, CATE, unlist(SCORE), return_df)
  names(return_list) = c("covariate_list", "matched_group", "match_quality", "matched_data")
  return(return_list)
}
