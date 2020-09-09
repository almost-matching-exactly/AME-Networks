require(Matrix)
require(igraph)
require(dplyr)
require(magrittr)
require(stringr)
require(Rcpp)
require(RcppArmadillo)
require(gmp)

my_fucking_ATE <- function(flout, Y, Z) {
  sapply(flout$MGs, function(x) {
    mean(Y[x][(Z[x] == 1)]) - mean(Y[x][(Z[x] == 0)])
  }) %>%
    mean()
}

my_combn <- function(x, m) {
  if (length(x) == 1) {
    return(list(x))
  }
  return(utils::combn(as.integer(x), m, simplify = FALSE))
}
# source('estimators.R')
# source('FLAME_bit.R')
# source('ATE.R')
# sourceCpp('subgraph_enumerate.cpp')

t1 <- graph(edges = c(1, 2), directed = FALSE)
t1 <- set_vertex_attr(t1, 'color', value = c(1, 0))
t2 <- graph(edges = c(1, 2), directed = FALSE)
t2 <- set_vertex_attr(t2, 'color', value = c(1, 1))

sq1 <- graph(edges = c(1, 2, 2, 3), directed = FALSE)
sq1 <- set_vertex_attr(sq1, 'color', value = c(1, 0, 0))
sq2 <- graph(edges = c(1, 2, 2, 3), directed = FALSE)
sq2 <- set_vertex_attr(sq2, 'color', value = c(1, 1, 0))
sq3 <- graph(edges = c(1, 2, 2, 3), directed = FALSE)
sq3 <- set_vertex_attr(sq3, 'color', value = c(1, 1, 1))

two_star1 <- graph(edges = c(), n = 2, directed = FALSE)
two_star1 <- set_vertex_attr(two_star1, 'color', value = c(1, 0))
two_star2 <- graph(edges = c(), n = 2, directed = FALSE)
two_star2 <- set_vertex_attr(two_star2, 'color', value = c(1, 1))

four_star1 <- graph(edges = c(), n = 4, directed = FALSE)
four_star1 <- set_vertex_attr(four_star1, 'color', value = c(1, 0, 0, 0))
four_star2 <- graph(edges = c(), n = 4, directed = FALSE)
four_star2 <- set_vertex_attr(four_star2, 'color', value = c(1, 1, 0, 0))
four_star3 <- graph(edges = c(), n = 4, directed = FALSE)
four_star3 <- set_vertex_attr(four_star3, 'color', value = c(1, 1, 1, 0))
four_star4 <- graph(edges = c(), n = 4, directed = FALSE)
four_star4 <- set_vertex_attr(four_star4, 'color', value = c(1, 1, 1, 1))

standardize <- function(x, type = '0-1', limits = c(0, 1)) {
  if (type == 'center') {
    return((x - mean(x, na.rm = T)) / (ifelse(sd(x) == 0, 1, sd(x))))
  } else if (type == '0-1') {
    if (length(unique(x)) == 1) {return(x)}
    return((limits[2] - limits[1]) / (max(x) - min(x)) * (x - max(x)) + limits[2])
  } else {
    stop("I don't know how to standardize this!")
  }
}

count_degree = function(G){
  dm = max(degree(G))
  X = matrix(0, length(V(G)), dm)
  for(i in V(G)){
    isg = induced_subgraph(G, neighbors(G, i))
    dgs = table(degree(isg))
    X[i, as.integer(names(dgs)) + 1] = dgs
  }
  X[, which(colSums(X) > 0)]
}

count_feature <- function(A, G, feature, Z) {
  # Takes in:
  #   The (treated) adjacency matrix A
  #   The (treated) graph G
  #   A string, feature, specifying which feature should be counted
  # Valid features are:
  #   'betweenness', 'closeness', 'degree', 'triangle',
  #   'kstar(n)' (which finds n-stars)
  # Returns:
  #   A vector of length n_units with the value of the feature for each unit
  n <- vcount(G)
  if (feature == 'betweenness') {
    return(betweenness(G, directed = FALSE, normalized = TRUE))
  }
  if (feature == 'closeness') {
    return(closeness(G, normalized = TRUE))
  }
  if (feature == 'degree') {
    return(colSums(A))
    return(colSums(A * Z))
  }
  else {
    feature_counts <- vector(mode = 'numeric', length = n)
  }
  for (i in 1:n) {
    isg <- induced_subgraph(G, neighbors(G, i))
    if (feature == 'triangle') {
      feature_counts[i] <- length(E(isg))
      # feature_counts[i] <- 
      #   count_subgraph_isomorphisms(t1, isg) + 
      #   count_subgraph_isomorphisms(t2, isg)
    } else if (feature == 'square') {
      feature_counts[i] <- 
        count_subgraph_isomorphisms(sq1, isg) + 
        count_subgraph_isomorphisms(sq2, isg) + 
        count_subgraph_isomorphisms(sq3, isg)
    } else if (feature == '2star') {
      feature_counts[i] <- 
        count_subgraph_isomorphisms(two_star1, isg) + 
        count_subgraph_isomorphisms(two_star2, isg)
    } else if (feature == '4star') {
      feature_counts[i] <- 
        count_subgraph_isomorphisms(four_star1, isg) + 
        count_subgraph_isomorphisms(four_star2, isg) +
        count_subgraph_isomorphisms(four_star3, isg) +
        count_subgraph_isomorphisms(four_star4, isg)
    } else if (feature == '3-degree-neighb') {
      degs <- which(degree(isg) >= 3)
      feature_counts[i] <- 0
      if (length(degs) == 0) {
        next
      }
      for (j in 1:length(degs)) {
        neighbs <- neighbors(isg, degs[j])
        if (sum(get.vertex.attribute(induced_subgraph(isg, neighbs), name = 'color')) > 0) {
          feature_counts[i] <- feature_counts[i] + 1
        }
      }
    }
    else {
      stop("I don't know what feature this is!")
    }
  }
  return(feature_counts)
}

threshold_list_subgraphs = function(V, k) {
  if (k == 'max') {
    k <- min(length(V), k)
  }
  else {
    k <- min(length(V), k)
  }
  sgs = c()
  for (j in 1:k) {
    if (j == 1 && length(V) == 1) {
      sgs = c(sgs, list(V)) # Can probably simplify
    } else if (j <= length(V)) {
      sgs = c(sgs, combn(V, j, simplify = FALSE))
    }
  }
  # for (v in combn(V, k, simplify = FALSE)){
  #   sgs = c(sgs, list.flatten(list_subgraphs(v)))
  # }
  sgs
}

threshold_all_neighborhood_subgraphs = function(G, k = 'max') {
  # For each vertex i, generates sgs[[i]] which is all possible combinations of
  # ≤ k vertices out of i's neighbors
  # k <- sort(ego_size(G), decreasing = TRUE)[2]
  sgs = list()
  for (i in V(G)) {
    if (i%%10 == 0) {
      print(i)
    }
    neighbs <- neighbors(G, i)
    if (length(neighbs) == 0)
      sgs[[i]] = numeric(0)
    else
      sgs[[i]] = threshold_list_subgraphs(neighbs, k)
  }
}

find_sg_type = function(G, sg, sg_types, coloring = coloring) {
  # Given a subgraph sg and sg_types, returns
  # i if sg is isomorphic to the i'th subgraph in sg_types
  # and NA otherwise
  if (length(sg_types) == 0)
    return(NA)
  sG = induced_subgraph(G, sg)
  coloring <- vertex_attr(sG, 'color')
  # There's at least one untreated individual
  # if (!is.null(coloring) & sum(coloring) < length(V(sG))) {
  #   return(-1)
  # }
  # speeds things up, if possible, when not coloring
  method <- ifelse(is.null(coloring), 'auto', 'vf2')
  for (i in 1:length(sg_types)) {
    if (is_isomorphic_to(sG, sg_types[[i]], method = method))
      return(i)
  }
  return(NA)
}

gen_features = function(G, sgs, sg_types = list(), coloring) {
  # sgs is a list, each entry of which is a subset of vertices in the
  # neighborhood of a single vertex in G
  feats = list()
  
  for (g in sgs) {
    type = find_sg_type(G, g, sg_types, coloring)
    
    # If we haven't yet seen this type of subgraph...
    if (is.na(type)) {
      # Introduce a new type and classify it as such
      type = length(sg_types) + 1
      sg_types[[type]] = induced_subgraph(G,g)
      feats[[type]] = 0
    }
    if (type == -1) { # wholly untreated subgraph
      next
    }
    if (type > length(feats) || is.null(feats[[type]]))
      feats[[type]] = 0
    
    # One more observation of this subgraph type
    feats[[type]] = feats[[type]] + 1
  }
  return(list(feats, sg_types))
}

gen_all_features = function(G, sgs, coloring) {
  sg_types = list()
  feats = list()
  for (i in 1:length(sgs)) { # For each unit
    res = gen_features(G, sgs[[i]], sg_types, coloring)
    sg_types = res[[2]]
    feats[[i]] = res[[1]]
  }
  return(list(feats = feats,
              types = sg_types))
}

gen_data = function(feats) {
  dta = NULL
  n_feats = max(sapply(feats, length))
  for (i in feats){
    lst = i
    if(length(lst) < n_feats)
      lst[[n_feats]] = 0
    vs = unlist(lapply(lst, function(x) ifelse(is.null(x), 0, x)))
    dta = rbind(dta, vs)
  }
  data.frame(dta, row.names = 1:length(feats))
}

new_sg_counter <- function(g, threshold) {
  n <- length(V(g))
  all_neighb_subgraphs <- vector(mode = 'list', length = n)
  for (i in 1:n) {
    tmp <- gen_connected_sgs(g, neighbors(g, i), sgs_so_far = c(), c(), threshold)
    all_neighb_subgraphs[[i]] <- tmp[2:length(tmp)]
  }
  return(all_neighb_subgraphs)
}

ATE_est_error <- function(Y, f, estimator_type, G, A, ATE, Z,
                          feature_counts, network_lik_weight, treat_prob, 
                          iterate_FLAME, threshold, COV=NULL) {
  # Takes in:
  #  The outcome vector Y
  #  The interference vector f
  #  The estimator type, one of:
  #    'true', 'first_eigenvector', 'all_eigenvectors', 'degree_dist', 'naive', 'FLAME'
  #  The graph G
  #  The (true) average treatment effect ATE
  # Returns absolute error between the estimated ATE and the true value, a scalar
  graph_dist <- NULL
  n <- length(Y)
  total_treatment_effect <- 0
  if (estimator_type == 'true') {
    graph_dist <- vector(mode = 'numeric', length = n)
    for (i in 1:n) {
      opposite_treatment <- which(Z != Z[i])
      mif <- abs(f - f[i])
      match <- opposite_treatment[which.min(mif[opposite_treatment])]
      Ymif <- Y[match]
      total_treatment_effect <- total_treatment_effect +
        (Y[i] - Ymif) * Z[i] + (Ymif - Y[i]) * (1 - Z[i])
      graph_dist[i] <- get_graph_dist(i, match, A)
    }
    graph_dist <- mean(graph_dist, na.rm = TRUE)
    error <- abs(total_treatment_effect / n - ATE)
    # error <- abs(true_est(Y, Z, f) - ATE)
  } else if (estimator_type == 'stratified') {
    error <- abs(strat_naive(A, Z, Y) - ATE)
  } else if (estimator_type == 'SANIA') {
    error <- abs(ind_SANIA(G, Z, Y, mean(Z)) - ATE)
  } else if (estimator_type == 'first_eigenvector') {
    error <- abs(first_eig(A, Z, Y) - ATE)
  } else if (estimator_type == 'all_eigenvectors') {
    eig <- eigen(A, symmetric = TRUE)
    total_treatment_effect <- 0
    ########## Might make more sense to weigh the evecs by their inverse evals, no...?
    weights <- diag(1 / (1:ncol(eig$vectors)))
    # graph_dist <- vector(mode = 'numeric', length = n)
    for (i in 1:n) {
      opposite_treatment <- which(Z != Z[i])
      mi <- sapply(1:n,
                   function(x) {
                     weighted_l2(eig$vectors[i, ], eig$vectors[x, ], weights)})
      match <- opposite_treatment[which.min(mi[opposite_treatment])]

      Ymie <- Y[match]

      total_treatment_effect <- total_treatment_effect +
        (Y[i] - Ymie) * Z[i] + (Ymie - Y[i]) * (1 - Z[i])

    #   # Graph distance
    #   ego_neighbors <- which(A[i, ] == 1)
    #   ego_neighborhood <- A[ego_neighbors, ego_neighbors]
    # 
    #   match_neighbors <- which(A[match, ] == 1)
    #   match_neighborhood <- A[match_neighbors, match_neighbors]
    # 
    #   if (length(ego_neighbors) == 0 & length(match_neighbors) == 0) {
    #     graph_dist[i] <- 0
    #     next
    #   }
    # 
    #   if (length(ego_neighbors) < length(match_neighbors)) {
    #     n_add <- length(match_neighbors) - length(ego_neighbors)
    #     ego_neighborhood <-
    #       bdiag(ego_neighborhood, matrix(0, nrow = n_add, ncol = n_add)) %>%
    #       matrix(nrow = nrow(.))
    #   } else if (length(ego_neighbors) > length(match_neighbors)) {
    #     n_add <- length(ego_neighbors) - length(match_neighbors)
    #     match_neighborhood <-
    #       bdiag(match_neighborhood, matrix(0, nrow = n_add, ncol = n_add)) %>%
    #       matrix(nrow = nrow(.))
    #   }
    #   if (any(dim(ego_neighborhood) == 0) |
    #       any(dim(match_neighborhood) == 0)) {browser()}
    #   graph_dist[i] <- minimal_frobenius(ego_neighborhood, match_neighborhood)
    #   # graph_dist[i] <- mean((ego_neighborhood - match_neighborhood) ^ 2)
    }
    # graph_dist <- mean(graph_dist, na.rm = TRUE)
    error <- abs(total_treatment_effect / n - ATE)
    # error <- abs(all_eigs(A, Z, Y) - ATE)
  } else if (estimator_type == 'FLAME') {
    
    ## Old method; avoid using if possible
    # all_subgraphs <- threshold_all_neighborhood_subgraphs(G, threshold)
    
    all_subgraphs <- get_neighb_subgraphs(A, 1:n, threshold)
    features_and_graphs <- gen_all_features(G, all_subgraphs)
    all_features <- features_and_graphs$feats
    sg_types <- features_and_graphs$types
    dta <- gen_data(all_features)
    
    if (!is.null(COV)) {
      dta <- cbind(dta, COV)  
    }
    dta <- data.frame(sapply(dta, factor), stringsAsFactors = T)
    dta$outcome <- Y
    dta$treated <- factor(Z)
    # browser()
    flame_out <- FLAME_bit(dta, dta, A = A,
                           network_lik_weight = network_lik_weight, iterate_FLAME = iterate_FLAME)
    
    # matches <- flame_out$MGs
    # total_graph_dist <- vector(mode = 'numeric', length = length(matches))
    # for (i in 1:length(matches)) {
    #   group <- matches[[i]]
    #   group_graph_dist <- vector(mode = 'numeric', length = length(group))
    #   for (j in 1:length(group)) {
    #     curr_unit <- group[j]
    #     opp_treatment <- group[Z[curr_unit] != Z[group]]
    #     graph_distances <- vector(mode = 'numeric', length = length(opp_treatment))
    #     for (k in 1:length(opp_treatment)) {
    #       graph_distances[k] <- get_graph_dist(curr_unit, opp_treatment[k], A)
    #     }
    #     group_graph_dist[j] <- mean(graph_distances, na.rm = TRUE)
    #   }
    #   total_graph_dist[i] <- mean(group_graph_dist, na.rm = TRUE)
    # }
    # graph_dist <- mean(total_graph_dist, na.rm = TRUE)
    #print(paste('Matched=', sum(flame_out$matched_data$matched)))
    #print(paste('ATE=', ATE(flame_out)))
    # error <- abs(ATE(flame_out) - ATE)
    error <- abs(my_fucking_ATE(flame_out, Y, Z) - ATE)
  } else if (estimator_type == 'naive') {
    error <- abs(mean(Y[Z == 1]) - mean(Y[Z == 0]) - ATE)
  } else {
    stop("I don't know what this estimator is!")
  }
  # print(sprintf('The graph distance of %s is %.4f', estimator_type, graph_dist))
  return(list(error = error,
              graph_dist = graph_dist))
}

minimal_frobenius <- function(A1, A2) {
  minimum <- Inf
  if (!is.matrix(A1) | !is.matrix(A2)) {
    size <- 1
  } 
  else {
    size <- nrow(A1)
  }
  
  if (size <= 1) {
    if (length(abs(c(A1 - A2))) == 0) {browser()}
    return(abs(c(A1 - A2)))
  }
  
  all_permutations <- gtools::permutations(size, size)
  for (i in 1:nrow(all_permutations)) {
    frob_norm <- 
      sqrt(sum((A1 - A2[all_permutations[i, ], ]) ^ 2, na.rm = TRUE) / size)
    if (frob_norm < minimum) {
      minimum <- frob_norm
    }
  }
  if (length(minimum) == 0) {browser()}
  return(minimum)
}

get_graph_dist <- function(unit1, unit2, A) {
  ego_neighbors <- which(A[unit1, ] == 1)
  ego_neighborhood <- A[ego_neighbors, ego_neighbors]
  
  match_neighbors <- which(A[unit2, ] == 1)
  match_neighborhood <- A[match_neighbors, match_neighbors]
  
  if (length(ego_neighbors) == 0 & length(match_neighbors) == 0) {
    return(0)
  }
  
  if (length(ego_neighbors) < length(match_neighbors)) {
    n_add <- length(match_neighbors) - length(ego_neighbors)
    ego_neighborhood <- 
      bdiag(ego_neighborhood, matrix(0, nrow = n_add, ncol = n_add)) %>% 
      matrix(nrow = nrow(.))
  } else if (length(ego_neighbors) > length(match_neighbors)) {
    n_add <- length(ego_neighbors) - length(match_neighbors)
    match_neighborhood <- 
      bdiag(match_neighborhood, matrix(0, nrow = n_add, ncol = n_add)) %>% 
      matrix(nrow = nrow(.))
  }
  if (any(dim(ego_neighborhood) == 0) | 
      any(dim(match_neighborhood) == 0)) {
    browser()
  }
  graph_dist <- minimal_frobenius(ego_neighborhood, match_neighborhood)
  return(graph_dist)
}

simulate_network_matching <- function(sim_type = 'ER',
                                      G0 = NULL,
                                      n_sims = 50,
                                      n_units = 50,
                                      n_clusters = 5,
                                      treat_prob = NULL,
                                      n_treated = NULL,
                                      interference_type = 'all',
                                      homoscedastic = TRUE,
                                      erdos_renyi_p = 0.07,
                                      standardization_type = 'center',
                                      interference_features,
                                      interference_parameters = NULL,
                                      matching_features = c(),
                                      estimators,
                                      network_lik_weight = 0,
                                      coloring = TRUE,
                                      n_blocks = 5,
                                      iterate_FLAME = FALSE,
                                      multiplicative = FALSE,
                                      covariate_weight = 0, 
                                      covariate_lvs = 3, 
                                      threshold) {
  
  if (network_lik_weight < 0) {
    stop('network_lik_weight should be positive, as it weighs the AIC')
  }
  
  stopifnot(length(interference_features) == length(interference_parameters))
  if ('degree_dist' %in% matching_features & length(matching_features) > 1) {
    stop('If matching on degree sequence, cannot match on other features as well')
  }
  
  # To store errors of the estimators
  abs_error <- matrix(, nrow = n_sims, ncol = length(estimators))
  graph_dist <- matrix(, nrow = n_sims, ncol = length(estimators))
  colnames(abs_error) <- estimators
  
  # Parameters used to generate outcome
  X = matrix(0, n_units, covariate_lvs)
  alpha = rep(0, n_units)
  for(i in 1:n_units) {
    X[i, sample(1:covariate_lvs, 1)] = 1
    alpha[i] = rnorm(1, which(X[i, ]==1) * covariate_weight, 1)
    # alpha[i] = rnorm(1, X[i, ] * c(-3, 4, 6) * covariate_weight, 1)
  }
  # if(covariate_weight==0) {
  #   X = NULL
  # }
  if (homoscedastic) {
    alpha <- rnorm(n_units, 0, 1)  
  }
  else {
    # alpha <- rnorm(n_units, 0, rep(runif(n_clusters, 0, sqrt(5)) ^ 2, each = n_clusters))
    alpha <- rnorm(n_units, 0, runif(n_units))
  }
  beta <- rnorm(n_units, 5, 1)
  ATE <- mean(beta)
  
  all_features <- base::union(interference_features, matching_features)
  n_features <- length(all_features)
  
  if (sim_type == 'SBM') {
    interference_parameters %<>%
      lapply(function(x) runif(n_blocks, x[1], x[2])) %>%
      unlist() %>%
      rep(each = n_units / n_blocks) %>%
      matrix(nrow = n_units)
    pmat = Matrix(runif(n_blocks ^ 2, 0, 0.1),
                  nrow = n_blocks,
                  ncol = n_blocks)
    pmat = forceSymmetric(pmat)
    bsizes = rep(n_units / n_blocks, n_blocks)
  }
  t0 = proc.time()[3]
  tot_time = 0
  for (sim in 1:n_sims) {
    t1 = proc.time()[3] - t0
    t0 = proc.time()[3]
    tot_time = tot_time + t1
    print(paste('Simulation', sim, 'of', n_sims, "last one took", round(t1, 1), 'seconds.', 
                'ETA:', round((tot_time/sim) * (n_sims - sim), 1), 'seconds.'))
    
    # Generate graph
    if (is.null(G0)) { # If no input graph
      if (sim_type == 'ER') {
        G <- erdos.renyi.game(n_units, erdos_renyi_p, directed = FALSE)
      } else if (sim_type == 'SBM') {
        G <- sbm.game(n_units, pmat, bsizes, directed = F, loops = F)
      } else if (sim_type == 'power law') {
        G <- barabasi.game(50, power = 2)
      } else if (sim_type == 'clustered') {
        # cluster_assignment <- 
        #   make_randomized_cluster(n_clusters = n_clusters, 
        #                           n_units = rep(n_units / n_clusters, n_clusters),
        #                           within_cluster_probs = rep(erdos_renyi_p, n_clusters),
        #                           n_treated = rep(n_units / n_clusters, n_clusters) / 2)
        # 
        # G <- cluster_assignment$graph
        pref_mat <- matrix(.005, nrow = n_clusters, ncol = n_clusters)
        diag(pref_mat) <- erdos_renyi_p
        G <- sbm.game(n_units, 
                      pref.matrix = pref_mat,
                      block.sizes = rep(n_units / n_clusters, n_clusters),
                      directed = FALSE,
                      loops = FALSE)
        Z <- vector(mode = 'numeric', length = sum(n_units))
        n_cum_units <- c(0, cumsum(rep(n_units / n_clusters, n_clusters)))
        for (i in 1:n_clusters) {
          Z[sample((n_cum_units[i] + 1):n_cum_units[i + 1], n_units / (2 * n_clusters))] <- 1
        }
        # Z <- cluster_assignment$treatment
      }
      else {
        stop("I don't recognize this type of simulation!")
      }
    } else {
      G <- G0
    }
    
    # Randomly assign treatment to n_treated units
    if ((is.null(n_treated) & is.null(treat_prob)) | 
        (!is.null(n_treated) & !is.null(treat_prob))){
      stop('You must specify one and only one of n_treated and treat_prob')
    }
    if (is.null(n_treated)) {
      # Assign treatment randomly with probability treat_prob
      Z <- sample(c(0, 1), size = n_units, replace = TRUE, prob = c(1 - treat_prob, treat_prob))
      n_treated <- sum(Z)
    }
    else {
      # Assign treatment randomly to n_treated units 
      if (sim_type != 'clustered') {
        Z <- rep(0, n_units)
        Z[sample(1:n_units, n_treated)] <- 1
      }
    }
    
    if (coloring) {
      G <- set_vertex_attr(G, 'color', value = Z)
    }
    # Get the graph adjacency matrix
    A <- get.adjacency(G, type= 'both', sparse = FALSE)
    feature_counts <- matrix(0, nrow = n_units, ncol = n_features)
    colnames(feature_counts) <- all_features
    
    if (interference_type == 'drop_mutual_untreated_edges') {
      AZ <- A
      untreated <- which(Z == 0)
      to_drop <- combn(untreated, 2)
      for (i in 1:ncol(to_drop)) {
        AZ[to_drop[1, i], to_drop[2, i]] <- 0
      }
      AZ %<>% forceSymmetric()
    } else if (interference_type == 'drop_untreated_edges') {
      untreated <- which(Z == 0)
      AZ <- A
      AZ[untreated, ] <- 0
      AZ[, untreated] <- 0
    } else {
      print('Assigning interference to the whole, untouched graph!')
      AZ <- A
    }
    GZ <- graph_from_adjacency_matrix(as.matrix(AZ), mode = 'undirected')
    for (i in 1:n_features) {
      feature_counts[, i] <- count_feature(AZ, GZ, all_features[i], Z)
    }
    std_feature_counts <- apply(feature_counts, 2, standardize, standardization_type)
    # check MARGIN!!
    if (sim_type == 'ER') {
      f <- sweep(std_feature_counts[, interference_features],
                 MARGIN = 2,
                 STATS = interference_parameters,
                 FUN = '*')
    } else {
      f <- std_feature_counts[, interference_features] * interference_parameters
    }
    f <- rowSums(f) 
    if (multiplicative) {
      tmp <- std_feature_counts[, interference_features]
      cts <- tmp[, 1] * tmp[, 2]
      f <- prod(interference_parameters) * cts
    }
    Y = alpha + beta * Z + f
    # resids <-
    #   lm(Y ~ X) %>%
    #   summary() %>%
    #   pluck('residuals')
    
    for (j in 1:length(estimators)) {
      # if (estimators[j] != 'FLAME') {
      #   Y <- resids
      # }
      estimation <- 
        ATE_est_error(Y, f, estimators[j], G, A, ATE, Z,
                      feature_counts, network_lik_weight, treat_prob, iterate_FLAME, threshold, COV = X)
      abs_error[sim, j] <- estimation$error
      # graph_dist[sim, j] <- estimation$graph_dist
    }
    print('Error:')
    print(abs_error[sim, ])
  }
  print(paste('Done, total time elapsed:', round(tot_time, 1)))
  mean_abs_error <- colMeans(abs_error)
  sd_abs_error <- apply(abs_error, 2, sd)
  
  # print(sprintf('For %d units, %s coloring....', n_units, c('without', 'with')[1 + coloring]))
  print(sprintf('With p = %.2f,', erdos_renyi_p))
  print(interference_parameters)
  print(interference_features)
  for (i in seq_along(estimators)) {
    'The mean absolute error of the %s estimator was %.2f' %>%
      sprintf(estimators[i], mean_abs_error[i]) %>%
      print()
  }
  return(abs_error)
  return(list(error = abs_error,
              graph_dist = graph_dist))
}