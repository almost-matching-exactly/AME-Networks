plot_results <- function(out_all, plot_title, sim_name, save = TRUE) {
  n_sims <- nrow(out_all[[1]]$results)
  n_estimators <- ncol(out_all[[1]]$results)
  
  ggdata = NULL
  for (i in 1:length(interference_params)){
    ggdata = rbind(ggdata, data.frame(as.data.frame(out_all[[i]]$results), setting=i))
  }
  
  ggdata = melt(ggdata, id.vars = 'setting')
  ggdata = ggdata[ggdata$variable!='true', ]
  ggdata$variable = factor(ggdata$variable,
                           levels=c('FLAME', 'first_eigenvector',
                                    'all_eigenvectors', 'naive',
                                    'stratified', 'SANIA'),
                           labels=c('FLAME-Networks', 'First\n Eigenvector',
                                    'All\n Eigenvectors', 'Naive',
                                    'Stratified', 'SANIA'))
  # ggdata$variable <- factor(ggdata$variable, levels = c('FLAME', 'true'), 
  #                           labels = c('FLAME-Networks', 'True\nInterference'))
  
  #this is for coloring
  leq = rep(NA, length(interference_params) * n_estimators)
  i = 1
  for(m in 1:n_estimators) {
    for(k in 1:length(interference_params)){
      means = tapply(ggdata[ggdata$setting==k, 'value'], 
                     ggdata[ggdata$setting==k, 'variable'], mean)
      leq[i] = (means[m] <= means['FLAME-Networks'])
      i = 1 + i
    }
  }
  repcol = rep(leq, each=n_sims)
  
  g <- 
  ggplot(ggdata[order(ggdata$variable), ], aes(x=variable, y=value)) + 
    geom_violin(aes(fill=repcol), draw_quantiles = 0.5 ) + 
    geom_hline(data = data.frame(y = tapply(ggdata[ggdata$variable=='FLAME-Networks', 'value'], 
                                            ggdata[ggdata$variable=='FLAME-Networks', 'setting'], mean), 
                                 setting=1:length(interference_params)), 
               aes(yintercept=y), linetype=2) + 
    scale_fill_brewer(palette = 'Set1') + 
    xlab('') + ylab('Absolute Error') + 
    ggtitle(paste('Simulation Results:', plot_title)) + 
    facet_grid(setting ~., scales = 'free_y') +
    theme_bw() + theme(legend.position='None', plot.title = element_text(hjust = 0.5), 
                       text = element_text(color='black', size=16)) 

  print(g)
  if (save) {
    ggsave(paste('Plots/', sim_name, '.png', sep=''), 
           width=10, height=8, units = 'in', device = 'png', dpi=300)
  }
  return(g)
}