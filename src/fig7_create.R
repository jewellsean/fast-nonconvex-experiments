## Figure 7 of Jewell et al. (2018) arXiv:1802.07380
## This figure compares the estimated spikes and calicum concentrations
## obtained through the unconstained and positive spike constrained 
## solution of eqns. (2) and (3), respectively, from the paper. 
## We use the spike finder dataset to illustrate the differences. 
## 
## In this script, we
## 1. Read and process spike finder data 
## 2. Estimate spikes and calcium concentrations corresponding to eqn. (2) and (3)
## 3. Plot results from (2) and (3) around a region of interest 
## 4. Save plot as pdf in specified location 
library(tidyverse)
library(latex2exp)
library(FastLZeroSpikeInference)

source("utils.R")

## data parameters
fps <- 100 # frames per second 
data_set <- 7 # spike finder dataset number 
ind <- 13 # cell number 
subset_is <- 1:20000 # subset of the data used to estimate the solution 
gam <- 0.98 # decay parameter 

## tuning parameters 
lambda_p2 <- 3.6
lambda_p3 <- 3.76

## file i/o paramters 
spikefinder_base_dir <- "/Users/jewellsean/Dropbox/research/fast_nonconvex_deconvolution_of_calcium_imaging_data/spike_finder_data/"
fig_save_dir <- "~/Desktop/figures/"

## plot parameters 
wd <- 0.1 # rectangle box width 

## 1. Read and process spike finder data
data_dir <- paste0(spikefinder_base_dir, "spikefinder.train/", data_set, ".train.")
calcium_dat <- readr::read_csv(paste0(data_dir, "calcium.csv")) %>% rename_all(col_renamer)
spike_dat <- readr::read_csv(paste0(data_dir, "spikes.csv")) %>% rename_all(col_renamer)
data <- subset_data(calcium_dat, spike_dat, subset_is)


## 2. Estimate spikes and calcium concentrations corresponding to eqn. (2) and (3)
fit_2 <- estimate_spikes(dat = data$calcium_d, gam = gam, lambda = lambda_p2, 
                         constraint = T, estimate_calcium = T)

fit_3 <- estimate_spikes(dat = data$calcium_d, gam = gam, lambda = lambda_p3, 
                         constraint = F, estimate_calcium = T)


## 3. Plot results from (2) and (3) around a region of interest and 
## 4. Save plot as pdf in specified location 
neg_spikes <- which(fit_3$estimated_calcium[fit_3$spikes] - gam * fit_3$estimated_calcium[fit_3$spikes - 1] < 0 )
neg_spikes_magnitude <- fit_3$estimated_calcium[fit_3$spikes[neg_spikes]] -gam * fit_3$estimated_calcium[fit_3$spikes[neg_spikes - 1]]
large_neg_spike_idx <- sort(neg_spikes_magnitude, index.return = T)$ix[1]

## plot around large negative spike
neg_spike_time <- fit_3$spikes[neg_spikes[large_neg_spike_idx]] * (1 / fps)
ts <- neg_spike_time + 2 * c(-1, 1)

## setup plot
times <- subset_is * (1 / fps)
true_spike <- which(data$spikes > 0) * (1 / fps)

pdf(paste0(fig_save_dir, "fig7.pdf"), height = 8, width = 16)
par(mfrow = c(1, 2), mar = c(5, 3, 4, 2) + 0.1)
plot_estimates(times, data$calcium_d, true_spike, fit_3$spikes, xlim = ts, xlab = "Time (s)")
lines(times, fit_3$estimated_calcium, col = "blue")
rect(neg_spike_time - wd, -1.25, neg_spike_time + wd, 10, border = "darkred", lwd = 2)


plot_estimates(times, data$calcium_d, true_spike, fit_2$spikes, xlim = ts, xlab = "Time (s)")
lines(times, fit_2$estimated_calcium, col = "blue")
rect(neg_spike_time - wd, -1.25, neg_spike_time + wd, 10, border = "darkred", lwd = 2)
dev.off()
