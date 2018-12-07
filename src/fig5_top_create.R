library(tidyverse)
require(magrittr)
require(yaml)

# global configs
args <- commandArgs(trailingOnly = TRUE)
config_file <- args[1]
configs <- yaml.load_file(config_file)

experiment_directory <- configs$output$local_directory

df <- NULL
experiment_files <- list.files(experiment_directory, pattern = ".csv")

if (is.null(experiment_files)) {
  stop(paste0('This script analyzes results from a large scale experiment. Results, generated from fig5_experiments.R,
              must be saved in experiment_directory: ', experiment_directory))
}

## collect resutls from all experiments
for (experiment_file in experiment_files) {
  file_i <- paste0(experiment_directory, "/", experiment_file)
  df <- rbind(df, read.csv(file_i))
}  

df <- df %>% mutate(l0_optim = ifelse(measure == "Cor", 1 - l0_optim, l0_optim), 
              l1_optim = ifelse(measure == "Cor", 1 - l1_optim, l1_optim))

df$measure <- factor(df$measure, levels = c("VR", "VP", "Cor"), labels = c("Van Rossum", "Victor-Purpura", "1 - Correlation"))

## Compare the solution of (1.3) and (3.18)
p <- df %>% filter(constrained == FALSE) %>% 
  ggplot(aes(x = l0_optim, y = l1_optim, color = as.factor(dataset))) + 
  geom_point() + 
  facet_wrap(~measure, scales = "free") + 
  geom_abline(intercept = 0, slope = 1) +
  ylab('Solution to (3.18)') +
  xlab('Solution to (1.3)') + 
  theme_bw() + 
  theme(legend.position="none", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(size = 1,  colour = "black")) 
filename <- paste0(configs$output$figure_directory, "fig5_compare_1-3_3-18.pdf")
ggsave(filename, p, height = 4, width = 8)


## Compare the solution of (1.2) and (3.18)
p <- df %>% filter(constrained == TRUE) %>% 
  ggplot(aes(x = l0_optim, y = l1_optim, color = as.factor(dataset))) + 
  geom_point() + 
  facet_wrap(~measure, scales = "free") + 
  geom_abline(intercept = 0, slope = 1) +
  ylab('Solution to (3.18)') +
  xlab('Solution to (1.2)') + 
  theme_bw() + 
  theme(legend.position="none", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(size = 1,  colour = "black")) 


filename <- paste0(configs$output$figure_directory, "fig5_compare_1-2_3-18.pdf")
ggsave(filename, p, height = 4, width = 8)

## Compare the solution to (1.2) and (1.3)

df_compare <- df %>% select(dataset, cell_i, measure, l0_optim, constrained) %>% spread(constrained, l0_optim) 
colnames(df_compare) <- c("dataset", "cell_i", "measure", "unconstrained", "constrained")

p <- df_compare %>% 
  ggplot(aes(x = unconstrained, y = constrained, color = as.factor(dataset))) + 
  geom_point() + 
  facet_wrap(~measure, scales = "free") + 
  geom_abline(intercept = 0, slope = 1) +
  ylab('Solution to (1.2)') +
  xlab('Solution to (1.3)') +
  theme_bw() + 
  theme(legend.position="none", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(size = 1,  colour = "black")) 

filename <- paste0(configs$output$figure_directory, "fig5_compare_1-2_1-3.pdf")
ggsave(filename, p, height = 4, width = 8)
