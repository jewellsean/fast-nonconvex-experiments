## Recreate Fig 3 and 4 of Jewell et al. (2018)
## First need to run experiments in fig3-4_experiments.R through script
## run_fig3-4_experiments_desktop.sh
library(dplyr)
library(ggplot2)
library(magrittr)
library(yaml)

source("utils/data_utils.R")

# CLI
args <- commandArgs(trailingOnly = TRUE)
config_file <- args[1]
configs <- yaml.load_file(config_file)

# get local directories
local_exp_directory <- configs$local_exp_directory
local_fig_directory <- configs$local_fig_directory


# read all experiment data
df <- read_csvs_in_directory(local_exp_directory)

## Figure 3: number of intervals in unconstrained and constrained problems 

out_summary <- 
  df %>% filter(alg %in% c("fast-uncon", "fast-poscon")) %>%
  mutate(alg = factor(alg, levels = c("fast-uncon", "fast-poscon"), 
                      labels = c("Algorithm 1", "Algorithm 2")), 
         poisMean = factor(poisMean, levels = c("0.001", "0.01", "0.1"), 
                           labels = c("Theta = 0.001", "Theta = 0.01", "Theta = 0.1"))) %>% 
  group_by(alg, n, poisMean) %>% 
  summarize(mean_n_intervals = mean(n_intervals), 
            se_n_intervals = sd(n_intervals) / sqrt(n())) %>% 
  mutate(mean_low_se = mean_n_intervals - se_n_intervals, 
         mean_high_se = mean_n_intervals + se_n_intervals)

p2 <- out_summary %>% ggplot(aes(x = n, y = mean_n_intervals, col = alg)) +
  geom_errorbar(aes(ymin = mean_low_se, ymax = mean_high_se)) + 
  facet_wrap(~alg + poisMean) + #, scales = "free") +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() + 
  xlab('Length of series (T)') + 
  ylab('Number of regions') + 
  theme(legend.position="none", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(size = 1,  colour = "black")) + 
  scale_color_manual(values=c("orange", "purple")) 

ggsave(p2, filename = paste0(local_fig_directory, "fig3.pdf"), height = 6, width = 9)



## Figure 4: timing comparisons 
dfSummary <- df %>% group_by(n, alg, poisMean) %>%
  mutate(timed = timed) %>%
  summarize(meanCompute = mean(timed), seCompute = sd(timed) / n(), N = n())

names <- c(
`0.001` = "Theta = 0.001",
`0.01` = "Theta = 0.01",
`0.1` = "Theta = 0.1")

p1 <- dfSummary %>%
ggplot(aes(x = n, y = meanCompute, color = alg)) +
    geom_line(size = 0.8) +
    geom_point() +
    xlab('Length of series (T)') +
    ylab('Computation time (s)') +
    scale_y_log10() +
    scale_x_log10() +
    theme_bw() +
    theme(legend.position="none", 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(size = 1,  colour = "black"),
      axis.text=element_text(size=14),
      axis.title=element_text(size=14),
      strip.text.x = element_text(size = 14)) +
    scale_color_manual(values=c("purple", "orange", "red","blue")) +
    facet_grid(. ~ poisMean, labeller = labeller(names))

ggsave(p1, filename = paste0(local_fig_directory, "fig4.pdf"), height = 3, width = 9)
