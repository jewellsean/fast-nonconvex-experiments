library(tidyverse)
require(magrittr)
require(yaml)
require(gridExtra)

# global configs
args <- commandArgs(trailingOnly = TRUE)
config_file <- args[1]
configs <- yaml.load_file(config_file)

experiment_directory <- configs$local_directory

df <- NULL
experiment_files_all <- list.files(experiment_directory, pattern = ".csv")

experiment_files <- c()

for (f in experiment_files_all) {
  if (length(grep(".+-constraint-TRUE-.+", f)) > 0) {
    experiment_files <- c(experiment_files, f)
  }
}


if (is.null(experiment_files)) {
  stop(paste0('This script analyzes results from a large scale experiment. Results, generated from fig5_experiments.R,
              must be saved in experiment_directory: ', experiment_directory))
}

## collect resutls from all experiments
for (experiment_file in experiment_files) {
  file_i <- paste0(experiment_directory, "/", experiment_file)
  
  info <- str_split(experiment_file, "-")[[1]]
  tmp <- read.csv(file_i)
  tmp$dataset <- as.numeric(info[2])
  tmp$cell_i <- as.numeric(info[4])
  tmp$measure <- info[6]
  
  df <- rbind(df, tmp)
}  



spike_avg <- function(spike_loc, spikes, window_size) {
    M <- length(spikes)
    low_ind <- max(1, spike_loc - window_size)
    high_ind <- min(spike_loc + window_size, M)
    real_spikes <- sum(spikes[low_ind:high_ind])
    return(real_spikes)
}

percent_it <- function(x) {rank(x)/length(x)}

window_size <- 5
(window_size * 2  + 1) * (1 / 100)


get_mags <- function(ex){
    nFittedSpikes <- sum(ex$spike_loc)
    mapping_onto_real_spikes <- numeric(nFittedSpikes)
    fittedSpikes <- which(ex$spike_loc > 0)

    for (i in 1:nFittedSpikes) {
        mapping_onto_real_spikes[i] <- spike_avg(fittedSpikes[i], ex$true_spikes, window_size)
    }

    out <- data.frame(mag_percentiles = percent_it(ex$spike_mag[fittedSpikes]),
    smooth_spikes = mapping_onto_real_spikes,
    dataset = ex$dataset[1],
    measure = ex$measure[1],
    cell_i = ex$cell_i[1])

    return(out)
}


unique_combinations <- df %>% filter(measure == "VR") %>% distinct(measure, dataset, cell_i)

out_df <- NULL

for (i in 1:dim(unique_combinations)[1]) {
    ex <- df %>% filter(measure == unique_combinations$measure[i],
    dataset == unique_combinations$dataset[i],
    cell_i == unique_combinations$cell_i[i])
    out_df <- rbind(out_df,
    get_mags(ex))
}


out_df %$%
plot(mag_percentiles, smooth_spikes, xlab = "Percentile of FastLZeroSpikes", ylab = "Actual # of Spikes over 0.1 s window", pch = 20, cex = 0.2)

p1 <- out_df %>%
ggplot(aes(x = mag_percentiles, y = smooth_spikes)) +
    geom_hex(aes(fill="#000000",alpha=log10(..count..)), fill="#0000ff") +
    scale_alpha_continuous("Log # of points", breaks=seq(0,10,1)) +
    geom_smooth(color = "black", lwd = 3, method = "loess") +
    xlab('Estimated spike magnitude (percentile)') +
    ylab('Actual number of spikes in 0.1 second window') +
    theme_bw() +
    theme(legend.position="none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(size = 1,  colour = "black"))


p2 <- out_df %>% mutate(dataset = as.factor(dataset)) %>%
    ggplot(aes(x = mag_percentiles, y = smooth_spikes, color = dataset)) +
    geom_smooth(aes(fill = dataset), lwd = 2, method = "loess") +
    xlab('Estimated spike magnitude (percentile)') +
    ylab('Actual number of spikes in 0.1 second window') +
    theme_bw() +
    theme(legend.position="none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(size = 1,  colour = "black"))



p <- grid.arrange(p1, p2, nrow = 1)

ggsave(p, filename = paste0(configs$figure_directory, "fig5_bottom.pdf"), height = 5, width = 10)



