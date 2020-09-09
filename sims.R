## The file from which to run simulations
## The function to call is simulate_network_matching, which takes in arguments:
#  sim_type = 'ER': either 'ER' for Erdos-Renyi graph simulation or 'SBM' for stochastic block model generation
#  n_sims = 50: the number of simulations to run
#  n_units = 50: the number of units in the network
#   Careful making this too large as it will result in a massive number of subgraphs / computation time
#  n_blocks = 5: the number of communities when performing SBM generation
#  n_treated = floor(n_units / 2): the number of treated units
#   for now, use treat_prob (below) instead 
#  treat_prob = 0.5: the probability that each unit is treated
#   Use this instead of n_treated because I've derived the SANIA estimator in this setting 
#  erdos_renyi_p = 0.07; the p parameter in ER graph simulation.
#   Careful making this too large as it will result in a massive number of subgraphs / computation time
#  standardization_type = 'center': how to standardize feature counts when producing interference
#   Either 'center' which standardizes to mean 0, sd 1, or '0-1' which standardizes to be in [0, 1]
#  estimators: the estimators to use. Any of:
#   'true': nearest neighbor on true interference
#   'naive': difference in means
#   'all_eigenvectors': weird eigenvector Mahalanobis-like nearest neighbor thing
#   'first_eigenvectors': the above but only with the largest eigenvector
#   'FLAME': our approach
#   'stratified': the Sussman Airoldi 2017 stratified naive estimator 
#   'SANIA': the Sussman Airoldi 2017 SANIA minimum variance linear unbiased estimator
#  coloring: whether to take treatment into account when testing for isomorphism
#   for now, let's just go ahead and leave this as FALSE, though results should be similar
#  network_lik_weight: how much to weigh the maximized network log likelihood when 
#   computing match quality. This should be positive.
#   Would recommend making this 0 until I find a better way to determine it; 
#   otherwise, will probably be useless at best and counterproductive at worst
#  interference_type = 'drop_mutual_untreated_edges': determines what features
#   to count when computing / assigning interference. 
#   ONLY SPECIFY THIS WHEN LOOKING AT MISSPECIFIED INTERFERENCE AS IN EXP 3O. 
#   In that case you may make it: 
#   'drop_mutual_untreated_edges': drops any edges between 2 control units
#  interference_features: the features to count when assigning interference. Any of:
#   'triangle' to count triangles
#   'kstar(n)' to count n-stars (sorry for ugly syntax; previously for ergm package compatibility)
#   'degree' to count degree
#   'betweenness' to compute vertex betweenness 
#   'closeness' to compute closeness centrality
#   'k-degree-neighb' to count number of neighbors with degree >= k (sorry for ugly syntax)
#  interference_parameters: the values by which to weight the above features.
#   For sim_type == 'ER':
#     A vector of numerics equal in length to length(interference_features)
#     e.g. interference_features = c('degree', 'kstar(2)') and interference_parameters = c(1, 2)
#       implies interference = 1 * degree count + 2 * 2-star count
#   For sim_type == 'SBM':
#     A list equal in length to length(interference_features)
#     Each entry of the list is a vector of length 2 with lower and upper bounds, respectively,
#       from which the weight for the respective interference feature is to be sampled uniformly
#     E.g. interference_features = c('degree', 'kstar(2)') and
#       interference_parameters = list(c(0, 1), c(2, 5)) implies
#         interference = gamma1 * degree count + gamma2 * 2-star count, where
#           gamma1 ~ U(0, 1) and gamma2 ~ U(2, 5)
#  iterate_flame = FALSE: 
#   a boolean for whether to perform FLAME (TRUE) or simply do 1-round of exact matching (FALSE)
#   the latter is much more efficient (obviously) and recommended for dense / large networks

# Setup -------------------------------------------------------------------

set.seed(42069)
# setwd('~/Dropbox/Duke/projects/learning_interference/Network-FLAME/../networks/Up to Date Code/')
setwd('~/Desktop/Up to Date Code/')
source('network_flame_sims.R')
source('plot.R')
require(ggplot2)
require(reshape2)
require(beepr)

additive_features <- list(c('degree', 'triangle', '2star', '4star',
                            '3-degree-neighb', 'betweenness', 'closeness'),
                          c('degree', 'triangle', '2star', '4star',
                            '3-degree-neighb', 'betweenness', 'closeness'),
                          c('degree', 'triangle', '2star', '4star',
                            '3-degree-neighb', 'betweenness', 'closeness'),
                          c('degree', 'triangle', '2star', '4star',
                            '3-degree-neighb', 'betweenness', 'closeness'))

additive_params <- list(c(0, 10, 0, 0, 0, 0, 0),
                        c(10, 10, 0, 0, 0, 0, 0),
                        c(0, 10, 1, 1, 1, 1, -1),
                        c(5, 1, 10, 1, 1, 1, -1))

multiplicative_features <- list(c('degree', 'triangle'), 
                              c('degree', 'betweenness'), 
                              c('triangle', 'betweenness'), 
                              c('triangle', '3-degree-neighb'))

multiplicative_params <- list(c(5, 1), c(5, 1), c(5, 1), c(5, 1))

all_estimators <- c('first_eigenvector', 'all_eigenvectors',
                    'naive',  'FLAME',  'stratified', 'SANIA')


# Section 3.1: Additive Interference --------------------------------------

sim_name <- '3.1_additive'

interference_features <- additive_features

interference_params <- additive_params

out_all <- vector(mode = 'list', length = length(interference_params))
settings = list(sim_type = 'ER', 
                n_sims = 50,
                erdos_renyi_p = 0.05,
                n_units = 50,
                n_treated = 25,
                estimators = all_estimators,
                threshold = 5)

for (i in 1:length(interference_params)) {
  print(paste('Setting', i, 'of', length(interference_params)))
  out_all[[i]] <- list(settings)
  out_all[[i]]$results <- do.call(simulate_network_matching, 
                                  c(list(
                                    interference_features = interference_features[[i]],
                                    interference_parameters = interference_params[[i]]), 
                                    settings))
}

save(out_all, file=paste('Simulation_data/', sim_name, '.RData', sep=""))
beep() 

plot_results(out_all, plot_title = 'Additive Interference', sim_name = sim_name)


# Section 3.2: Covariate Adjustment ---------------------------------------

sim_name <- '3.2_cov_adjustment_weight_'    
cov_weights = c(5, 10, 15, 20, 25)
n_settings = 5
out_all <- vector(mode = 'list', length = n_settings)
settings = list(sim_type = 'ER', 
                n_sims = 50,
                n_units = 50,
                n_treated = 25,
                erdos_renyi_p = 0.05,
                estimators = all_estimators,
                interference_features = c('triangle', 'degree', 'betweenness'),
                interference_parameters = c(3, 1, 1),
                threshold = 5)
for (i in 1:n_settings) {
  print(paste('Setting', i, 'of', n_settings))
  out_all[[i]] <- list(settings)
  out_all[[i]]$results <- do.call(simulate_network_matching, 
                                  c(list(covariate_weight=cov_weights[i]), settings))
  save(out_all, file=paste('Simulation_data/', sim_name, i, '.RData', sep=""))
}
beep()

#### Plotting code
n_estimators <- ncol(out_all[[1]]$results)
ggdata = data.frame(matrix(NA, n_estimators * n_settings, 5))
names(ggdata) = c('method', 'setting', 'mean', 'lo', 'hi')
k = 1
for (i in 1:n_settings){
  for(j in 1:n_estimators){
    m = out_all[[i]]$results[, j]
    ggdata[k, 'method'] = colnames(out_all[[i]]$results)[j]
    ggdata[k, 'setting'] = cov_weights[i]
    ggdata[k, 'mean'] = mean(m)
    ggdata[k, 'lo'] = sort(m)[0.25 * length(m)]
    ggdata[k, 'hi'] = sort(m)[0.75 * length(m)]
    k = k + 1
  }
}

ggdata$method = factor(ggdata$method, levels = c('FLAME', 'first_eigenvector', 
                                                 'all_eigenvectors', 'naive', 'SANIA', 'stratified', 'true'), 
                       labels= c('FLAME-Networks', 'First\n Eigenvector', 
                                 'All Eigenvectors', 'Naive',  'SANIA', 'Stratified', 'True'))

ggplot(ggdata[!ggdata$method %in% c('First\n Eigenvector'), ], 
       aes(x=setting, y=mean, ymin=lo, ymax=hi, color=method, linetype=method)) + 
  geom_line(size=1) + geom_point() + 
  scale_color_brewer(palette = 'Set1') + 
  labs(color='Method', linetype='Method') + 
  xlab(TeX('Importance of Covariate ($\\beta$)')) + ylab('Mean Absolute Error') + 
  ggtitle('Simulation Results: Covariate Adjustment') + 
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), 
                     text = element_text(color='black', size=16),
                     legend.position = c(0.15, 0.8), legend.background = element_rect(color='black')) 

ggsave(paste('Plots/', sim_name, '.png', sep=''), width=10, height=7, units = 'in', device = 'png', dpi=300)


# Section 3.3: Misspecified Interference ----------------------------------

sim_name <- '3.3_misspecified'
interference_features <- c('degree', 'triangle')
interference_params <- list(c(5, 0), c(4, 1), c(3, 2), c(2, 3), c(1, 4), c(0, 5))

out_all <- vector(mode = 'list', length = length(interference_params))
settings = list(sim_type = 'ER', 
                n_sims = 50,
                n_units = 75,
                n_treated = 37,
                erdos_renyi_p = 0.05,
                interference_type = 'drop_mutual_untreated_edges',
                estimators = c('FLAME', 'stratified', 'SANIA'),
                interference_features = interference_features,
                threshold = 5)
for (i in 1:length(interference_params)) {
  print(paste('Setting', i, 'of', length(interference_params)))
  out_all[[i]] <- list(settings)
  out_all[[i]]$results <- do.call(simulate_network_matching, 
                                  c(list(interference_parameters=interference_params[[i]]), 
                                    settings))
  tmp <- out_all[[i]]
  save(tmp, file=paste('Simulation_data/', sim_name, '_tmp_', i, '.RData', sep=""))
}
beep()
save(out_all, file=paste('Simulation_data/', sim_name, '.RData', sep=""))

ggdata = NULL
for (i in 1:length(interference_params)){
  ggdata = rbind(ggdata, data.frame(as.data.frame(out_all[[i]]$results), setting=i))
}

ggdata = data.frame(mean = c(tapply(ggdata$FLAME, ggdata$setting, median),
                             tapply(ggdata$stratified, ggdata$setting, median),
                             tapply(ggdata$SANIA, ggdata$setting, median)),
                    lo = c(tapply(ggdata$FLAME, ggdata$setting,
                                  function(x) sort(x)[length(x) * 0.25]),
                           tapply(ggdata$stratified, ggdata$setting,
                                  function(x) sort(x)[length(x) * 0.25]),
                           tapply(ggdata$SANIA, ggdata$setting,
                                  function(x) sort(x)[length(x) * 0.25])),
                    hi = c(tapply(ggdata$FLAME, ggdata$setting,
                                  function(x) sort(x)[length(x) * 0.75]),
                           tapply(ggdata$stratified, ggdata$setting,
                                  function(x) sort(x)[length(x) * 0.75]),
                           tapply(ggdata$SANIA, ggdata$setting,
                                  function(x) sort(x)[length(x) * 0.75])),
                    setting = rep(0:5, 3), method=rep(c('FLAME-Networks', 'Stratified', 'SANIA'), each=6))

ggplot(ggdata, aes(x=setting, y=mean, color=method, fill=method, ymin=lo, ymax=hi)) + 
  geom_ribbon(alpha=0.2, color=NA) + geom_line(size=1) +  geom_point() + 
  xlab(TeX('$\\gamma$')) + ylab('Mean Absolute Error') + 
  scale_color_brewer(palette='Set1') + scale_fill_brewer(palette='Set1') + 
  theme_bw() + theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5), 
                     text = element_text(color='black', size=20), 
                     legend.position = c(0.82,0.89), 
                     legend.background = element_rect(colour = 'black')) 

beep()
ggsave(paste('Plots/3.3_misspecified', '.png', sep=''), width=10, height=7.5, units = 'in', device = 'png', dpi=300)

# Appendix G: Multiplicative Interference ---------------------------------

sim_name <- 'G_multiplicative'

interference_features <- multiplicative_features
interference_params <- multiplicative_params

out_all <- vector(mode = 'list', length = length(interference_params))
settings = list(sim_type = 'ER', 
                n_sims = 50,
                n_units = 50,
                n_treated = 25,
                erdos_renyi_p = 0.05,
                estimators = all_estimators,
                multiplicative = TRUE, 
                threshold = 5)
for (i in 1:length(interference_params)) {
  print(paste('Setting', i, 'of', length(interference_params)))
  out_all[[i]] <- list(settings)
  out_all[[i]]$results <- do.call(simulate_network_matching, 
                                  c(list(interference_features = interference_features[[i]], 
                                         interference_parameters=interference_params[[i]]), 
                                    settings))
}

save(out_all, file=paste('Simulation_data/', sim_name, '.RData', sep=""))
beep() 

plot_results(out_all, plot_title = 'Multiplicative Interference', sim_name = sim_name)

# Appendix H: Graph Cluster Randomization ---------------------------------
sim_name <- 'H_graph_cluster'

out_all <- vector(mode = 'list', length = length(interference_params))
settings = list(sim_type = 'clustered', 
                n_sims = 50,
                erdos_renyi_p = 0.05,
                n_units = 50,
                n_clusters = 5,
                n_treated = 25,
                estimators = c('first_eigenvector', 'all_eigenvectors',
                               'naive',  'FLAME',  'stratified', 'SANIA'),
                threshold = 5)
for (i in 1:length(interference_params)) {
  print(paste('Setting', i, 'of', length(interference_params)))
  out_all[[i]] <- list(settings)
  out_all[[i]]$results <- do.call(simulate_network_matching, 
                                  c(list(
                                    interference_features = interference_features[[i]],
                                    interference_parameters = interference_params[[i]]), 
                                    settings))
}

save(out_all, file=paste('Simulation_data/', sim_name, '.RData', sep=""))
beep() 

plot_results(out_all, plot_title = 'Graph Cluster Randomization', sim_name = sim_name)

# Appendix I: Real Network (AddHealth) ------------------------------------

sim_name <- 'I_real_network'

interference_features <- additive_features
interference_params <- additive_params

library(amen)
G <- 
  graph_from_adjacency_matrix(addhealthc3$Y, mode = 'undirected') %>%
  simplify()

out_all <- vector(mode = 'list', length = length(interference_params))
settings = list(G0 = G,
                n_sims = 50,
                n_units = 32,
                n_clusters = 5,
                n_treated = 16,
                estimators = all_estimators,
                threshold = 3)
for (i in 1:length(interference_params)) {
  print(paste('Setting', i, 'of', length(interference_params)))
  out_all[[i]] <- list(settings)
  out_all[[i]]$results <- do.call(simulate_network_matching, 
                                  c(list(
                                    interference_features = interference_features[[i]],
                                    interference_parameters = interference_params[[i]]), 
                                    settings))
  tmp <- out_all[[i]]
  save(tmp, file=paste('Simulation_data/', sim_name, '_tmp_', i, '.RData', sep=""))
}

save(out_all, file=paste('Simulation_data/', sim_name, '.RData', sep=""))
beep() 

plot_results(out_all, plot_title = 'AddHealth Network', sim_name = sim_name)

# Appendix J: Heteroscedastic Errors --------------------------------------


sim_name <- 'J_heteroscedastic'

interference_features <- additive_features
interference_params <- additive_params

out_all <- vector(mode = 'list', length = length(interference_params))
settings = list(sim_type = 'ER', 
                n_sims = 50,
                erdos_renyi_p = 0.07,
                n_units = 50,
                estimators = all_estimators,
                n_treated = 25,
                homoscedastic = FALSE,
                threshold = 5)
for (i in 1:length(interference_params)) {
  print(paste('Setting', i, 'of', length(interference_params)))
  out_all[[i]] <- list(settings)
  out_all[[i]]$results <- do.call(simulate_network_matching, 
                                  c(list(
                                    interference_features = interference_features[[i]],
                                    interference_parameters = interference_params[[i]]), 
                                    settings))
}

save(out_all, file=paste('Simulation_data/', sim_name, '.RData', sep=""))
beep() 

plot_results(out_all, plot_title = 'Heteroskedasticity', sim_name = sim_name)

# Appendix K: True Interference -------------------------------------------

sim_name <- 'K_true_interference'

interference_features <- additive_features
interference_params <- additive_params

out_all <- vector(mode = 'list', length = length(interference_params))
settings = list(sim_type = 'ER', 
                n_sims = 50,
                erdos_renyi_p = 0.06,
                n_units = 50,
                n_treated = 25,
                estimators = c('true', 'FLAME'),
                threshold = 5)
for (i in 1:length(interference_params)) {
  print(paste('Setting', i, 'of', length(interference_params)))
  out_all[[i]] <- list(settings)
  out_all[[i]]$results <- do.call(simulate_network_matching, 
                                  c(list(
                                    interference_features = interference_features[[i]],
                                    interference_parameters = interference_params[[i]]), 
                                    settings))
  tmp <- out_all[[i]]
  save(tmp, file=paste('Simulation_data/', sim_name, '_tmp_', i, '.RData', sep=""))
}

save(out_all, file=paste('Simulation_data/', sim_name, '.RData', sep=""))
beep() 

plot_results(out_all, plot_title = 'FLAME-Networks vs. Matching on True Interference', 
             sim_name = sim_name)