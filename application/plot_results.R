require(igraph)
source('network_flame_sims.R')
source('FLAME_bit.R')
source('ATE.R')
library(magrittr)
library(dplyr)
require(readxl)
require(tidyr)
library(ggplot2)
library(stringr)
library(reshape2)

set.seed(42069)

######### IMPORTANT: CHANGE TO YOUR DIRECTORY BEFORE RUNNING ############
setwd('/Users/mm567/Dropbox/Duke/projects/learning_interference/networks/Up to Date Code/application/')

FLAME_net = matrix(NA, 77, 3)
for (v_path in list.files('application_output/matches/')){
  pars = str_extract_all(v_path, "[0-9]+")[[1]]
  vq = readRDS(paste('application_output/matches/', v_path, sep=""))
  effs = unlist(sapply(vq$matched_group, "[", 'effect'))
  sizes = unlist(sapply(vq$matched_group, "[", 'size'))
  te = sum(effs * sizes)/sum(sizes)
  FLAME_net[as.integer(pars[1]), as.integer(pars[2])] = te
}

naive = matrix(NA, 77, 3)
for (v_path in list.files('Data/characteristics_1/')){
  vn = as.integer(str_extract(v_path, "[0-9]+"))
  vq1 = read.csv(paste('Data/characteristics_1/', v_path, sep=""))
  vq2 = read.csv(paste('Data/characteristics_2/', v_path, sep=""))
  vq3 = read.csv(paste('Data/characteristics_3/', v_path, sep=""))
  naive[vn, 1] = mean(vq1$Y[vq1$Z==1]) - mean(vq1$Y[vq1$Z==0])
  naive[vn, 2] = mean(vq2$Y[vq2$Z==1]) - mean(vq2$Y[vq2$Z==0])
  naive[vn, 3] = mean(vq3$Y[vq3$Z==1]) - mean(vq3$Y[vq3$Z==0])
}

diffs = naive - FLAME_net

FLAME_net = FLAME_net %>% data.frame() %>% melt() %>% drop_na
diffs = diffs %>% data.frame() %>% melt() %>% drop_na
naive = naive %>% data.frame() %>% melt() %>% drop_na

levels(diffs$variable) <- c('1', '2', '3')
levels(FLAME_net$variable) <- c('1', '2', '3')
levels(naive$variable) <- c('1', '2', '3')

boxplot_width <- 0.9
pos <- position_dodge2(preserve = 'total')
p_FLAME <- 
  ggplot(FLAME_net, aes(y = value, x = variable, color = variable, fill = variable)) +
  geom_boxplot(alpha = 0.4, width = boxplot_width, position = pos) +
  scale_color_brewer(palette = 'Set1') +
  scale_fill_brewer(palette = 'Set1') +
  labs(x = 'Question Number', y = 'FLAME-Network ADE Estimate', 
       fill = 'Question Number', color = 'Question Number') +
  theme_bw() +
  theme(legend.position = 'none',
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylim(c(-0.85, 0.3))

p_naive <- 
  ggplot(naive, aes(y = value, x = variable, color = variable, fill = variable)) +
  geom_boxplot(alpha = 0.4, width = boxplot_width, position = pos) +
  scale_color_brewer(palette = 'Set1') +
  scale_fill_brewer(palette = 'Set1') +
  labs(x = 'Question Number', y = 'Naive ADE Estimate', fill = 'Question Number', color = 'Question Number') +
  theme_bw() +
  theme(legend.position = 'none',
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylim(c(-0.7, 0.3))

p_diffs <- 
  ggplot(diffs, aes(y = value, x = variable, color = variable, fill = variable)) +
  geom_boxplot(alpha = 0.4, width = boxplot_width, position = pos) +
  scale_color_brewer(palette = 'Set1') +
  scale_fill_brewer(palette = 'Set1') +
  labs(x = 'Question Number', y = 'Difference in ADE Estimates:\n Naive - FLAME-Network', fill = 'Question Number', color = 'Question Number') +
  theme_bw() +
  theme(legend.position = 'none',
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylim(c(-.15, .2))

ggsave('application_output/plots/naive.png', plot = p_naive, 
       width=5, height=7, units = 'in', device = 'png', dpi=300)
ggsave('application_output/plots/FLAME.png', plot = p_FLAME, 
       width=5, height=7, units = 'in', device = 'png', dpi=300)
ggsave('application_output/plots/diff.png', plot = p_diffs, 
       width=5, height=7, units = 'in', device = 'png', dpi=300)
